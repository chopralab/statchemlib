#include "KBDynamics.hpp"

#include <thread>
#include "statchem/helper/logger.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "statchem/fileio/fileout.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/modeler/forcefield.hpp"
#include "statchem/modeler/modeler.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
#include "statchem/score/kbff.hpp"

#include "programs/common.hpp"

using namespace statchem_prog;
using namespace std;

template <>
ProgramInfo statchem_prog::program_information<KBDynamics>() {
    return ProgramInfo("kb_dynamics")
        .description(
            "Run dynamics on a complex using a knowledge-based potential.");
}

KBDynamics::KBDynamics() {}

bool KBDynamics::process_options(int argc, char* argv[]) {
    auto starting_inputs = common_starting_inputs();

    auto score_options = scoring_options();
    score_options.add_options()(
        "step_size", po::value<double>(&__step_size)->default_value(0.01),
        "Step size for fitted bsplines.");

    auto ff_min = forcefield_options();
    ff_min.add_options()("distance_cutoff",
                         po::value<double>(&__dist_cut)->default_value(6.0),
                         "Distance cutoff for intermolecular forces.")(
        "scale", po::value<double>(&__scale)->default_value(1.0),
        "Scale factor for the knowledge-based force.");

    auto dynamics_options = po::options_description("Dynamics Opetions");
    dynamics_options.add_options()(
        "temperature",
        po::value<double>(&__temperature)->default_value(300.0, "300"),
        "Temperature to run the dynamic simulation at.")(
        "friction", po::value<double>(&__friction)->default_value(91.0, "91.0"),
        "Friction/Collision frequency for a dynamics simulation in 1/ps")(
        "integrator",
        po::value<std::string>(&__integrator)->default_value("verlet"),
        "Integrator to use: 'verlet', 'langevin', or 'brownian'")(
        "dynamic_steps",
        po::value<int>(&__dynamics_steps)->default_value(1000))(
        "dynamic_step_size",
        po::value<double>(&__dynamics_step_size)->default_value(2.0, "2.0"),
        "Step size (in fempto seconds)");

    auto openmm = openmm_options();

    po::options_description cmdln_options;
    cmdln_options.add(starting_inputs);
    cmdln_options.add(score_options);
    cmdln_options.add(ff_min);
    cmdln_options.add(dynamics_options);
    cmdln_options.add(openmm);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdln_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        __help_text << "This program performs energy molecular dynamics using "
                       "bonded terms from Amber forcefield and nonbonded terms "
                       "from a generalized statistical potential.\n\n";
        __help_text << starting_inputs_help() << "\n";
        __help_text << cmdln_options << std::endl;
        return false;
    }

    __constant_receptor =
        process_starting_inputs(vm, __receptor_mols, __ligand_mols);

    process_scoring_options(vm, __ref, __comp, __func, __cutoff, __dist);

    if (__dist_cut > __cutoff) {
        throw std::out_of_range(
            "The --distance_cutoff must be <= the scoring function cutoff.");
    }

    process_forcefield_options(vm, __ffield, __mini_tol, __iter_max);
    process_openmm_options(vm, __platform, __precision, __accelerators, __checkpoint);

    return true;
}

int KBDynamics::run() {

    if (__receptor_mols.get_idatm_types().size() == 1) {
        __receptor_mols.compute_idatm_type()
            .compute_hydrogen()
            .compute_bond_order()
            .compute_bond_gaff_type()
            .refine_idatm_type()
            .erase_hydrogen()
            .compute_hydrogen()
            .compute_ring_type()
            .compute_gaff_type()
            .erase_hydrogen();
    }

    if (__ligand_mols.get_idatm_types().size() == 1) {
        __ligand_mols.compute_idatm_type()
            .compute_hydrogen()
            .compute_bond_order()
            .compute_bond_gaff_type()
            .refine_idatm_type()
            .erase_hydrogen()
            .compute_hydrogen()
            .compute_ring_type()
            .compute_gaff_type()
            .erase_hydrogen();
    }

    __score = std::unique_ptr<statchem::score::KBFF>(new statchem::score::KBFF(
        __ref, __comp, __func, __cutoff, __step_size));

    __score
        ->define_composition(__receptor_mols.get_idatm_types(),
                             __ligand_mols.get_idatm_types())
        .process_distributions(__dist)
        .compile_scoring_function();
    __score->compile_objective_function();

    __ffield.add_kb_forcefield(*__score, __dist_cut);

    statchem::OMMIface::SystemTopology::loadPlugins();

    /* 
        Need to hold track of the i and j variables in the for loops below.
        OpenMM state does not do that for us, so we use a separte file to store the contents.
    */
    int x, y;
    if(__checkpoint != "") {
        std::string candock_checkpoint = __checkpoint + "_candock";
        std::ifstream file(candock_checkpoint, ios::in);
        if(file.is_open()){
            file >> x;
            file >> y;
            file.close();
        }
        else {
            cerr << "Error could not open file: " << candock_checkpoint << endl;
            exit(1);
        }
    }
    
    for (size_t i = 0; i < __ligand_mols.size(); ++i) {

        if(__checkpoint != "")
            i = x;
    
        statchem::molib::Molecule& protein =
            __constant_receptor ? __receptor_mols[0] : __receptor_mols[i];
        statchem::molib::Molecule& ligand = __ligand_mols[i];

        statchem::molib::Atom::Grid gridrec(protein.get_atoms());
        protein.prepare_for_mm(__ffield, gridrec);

        __ffield.insert_topology(protein);
        __ffield.insert_topology(ligand);

        statchem::OMMIface::Modeler modeler(
            __ffield, "kb", __scale, __mini_tol, __iter_max, false,
            __dynamics_step_size, __temperature, __friction, __cutoff);

        modeler.set_num_steps_to_run(__dynamics_steps);

        modeler.add_topology(protein.get_atoms());
        modeler.add_topology(ligand.get_atoms());

        modeler.init_openmm(
            __platform, __precision, __accelerators,
            statchem::OMMIface::SystemTopology::integrator_type::verlet);

        modeler.add_crds(protein.get_atoms(), protein.get_crds());
        modeler.add_crds(ligand.get_atoms(), ligand.get_crds());

        modeler.unmask(ligand.get_atoms());
        modeler.unmask(protein.get_atoms());

        modeler.init_openmm_positions();

        modeler.minimize_state();

        modeler.__system_topology.set_temperature();
        modeler.__system_topology.set_box_vector();


        
        size_t j = 0;
        // Load checkpoint
        if(__checkpoint != "") {
            std::cerr << "Loading from checkpoint\n";
            modeler.__system_topology.load_checkpoint(__checkpoint);
            __checkpoint = "";

            j = y;
        }
        
        for (; j < 100; j++) {
            std::cerr << "\nligand_mols " << i << std::endl;
            std::cerr << "dynamics iteration " << j << std::endl;
            modeler.dynamics();
            modeler.__system_topology.save_checkpoint_candock(i, j);

            // init with minimized coordinates
            statchem::molib::Molecule minimized_receptor(
                protein, modeler.get_state(protein.get_atoms()));
            statchem::molib::Molecule minimized_ligand(
                ligand, modeler.get_state(ligand.get_atoms()));

            minimized_receptor.undo_mm_specific();

            statchem::fileio::print_complex_pdb(std::cout, minimized_ligand,
                                                minimized_receptor, 0.000);

            minimized_receptor.prepare_for_mm(__ffield, gridrec);
        }
        __ffield.erase_topology(ligand);
    }

    return 0;
}
