#include "PhysMinimize.hpp"

#include <thread>

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

using namespace statchem_prog;

template <>
ProgramInfo statchem_prog::program_information<PhysMinimize>() {
    return ProgramInfo("phys_min")
        .description("Minimize a complex using a physics-based potential.");
}

PhysMinimize::PhysMinimize() {}

int PhysMinimize::run(int argc, char* argv[]) {
    po::options_description starting_inputs("Starting input files");
    starting_inputs.add_options()("help,h", "Show this help menu.")(
        "receptor,r", po::value<std::string>()->default_value("receptor.pdb"),
        "Receptor filename. Will be ignored if the complex option is given")(
        "ligand,l", po::value<std::string>()->default_value("ligand.mol2"),
        "Ligand filename. Will be ignored if the complex option is given")(
        "complex,c", po::value<std::string>(),
        "File with both the ligand and receptor present. Overrides above "
        "options");

    po::options_description ff_min("Forcefield and Minimization Options");
    ff_min.add_options()("amber_xml", po::value<std::string>()->default_value(
                                          "data/amber10.xml"),
                         "Receptor XML parameters (and topology) input file")(
        "water_xml", po::value<std::string>()->default_value("data/tip3p.xml"),
        "Water XML parameters (and topology) input file")(
        "gaff_dat", po::value<std::string>()->default_value("data/gaff.dat"),
        "Gaff DAT forcefield input file")(
        "mini_tol", po::value<double>()->default_value(0.0001),
        "Minimization tolerance")(
        "max_iter", po::value<int>()->default_value(10),
        "Maximum iterations for minimization during linking");

    po::options_description openmm("OpenMM options");
    openmm.add_options()(
        "platform", po::value<std::string>()->default_value("CPU"),
        "Platform to run KBForce on. Options are CPU, GPU, and OpenCL.")(
        "precision", po::value<std::string>()->default_value("double"),
        "Precision to run KBForce on. Options are single, mixed, double. "
        "Only works using CUDA or OpenCL platform")(
        "accelerators", po::value<std::string>()->default_value("double"),
        "Precision to run KBForce on. Options are single, mixed, double. "
        "Only works using CUDA or OpenCL platform");

    po::options_description cmdln_options;
    cmdln_options.add(starting_inputs);
    cmdln_options.add(ff_min);
    cmdln_options.add(openmm);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdln_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout
            << "This program is used to minimize a given complex\n"
               "using a physics-based forcefield. To use it, specify a "
               "'receptor' and a ligand'\n"
               "file. Alternatively, you can specify a complex file and\n"
               "the residue name of interest. You may also specify files\n"
               "containing multiple complexes. The following cases are valid\n"
            << std::endl;
        std::cout << "(1) A single complex PDB file with multiple complexes "
                     "separted by MODEL and ENDMDL"
                  << std::endl;
        std::cout << "(2) A receptor file containing a single biomolecule and "
                     "a ligand file containing mutliple small-molecules"
                  << std::endl;
        std::cout << "(3) A receptor file containing multiple biomolecules and "
                     "a ligand file containing the same number of "
                     "small-molecules"
                  << std::endl;
        std::cout << std::endl << "Valid options are: " << std::endl;
        std::cout << cmdln_options << std::endl;

        return 0;
    }

    statchem::molib::Molecules receptor_mols;
    statchem::molib::Molecules ligand_mols;

    bool constant_receptor = false;
    if (vm.count("complex")) {
        auto complex = vm["complex"].as<std::string>();

        statchem::parser::FileParser drpdb(
            complex, statchem::parser::pdb_read_options::protein_poses_only |
                         statchem::parser::pdb_read_options::all_models);

        drpdb.parse_molecule(receptor_mols);

        statchem::parser::FileParser dlpdb(
            complex, statchem::parser::pdb_read_options::docked_poses_only |
                         statchem::parser::pdb_read_options::skip_atom |
                         statchem::parser::pdb_read_options::all_models);
        dlpdb.parse_molecule(ligand_mols);

    } else {
        auto receptor = vm["receptor"].as<std::string>();
        statchem::parser::FileParser rpdb(
            receptor, statchem::parser::pdb_read_options::all_models);
        rpdb.parse_molecule(receptor_mols);

        auto ligand = vm["ligand"].as<std::string>();
        statchem::parser::FileParser lpdb(
            ligand, statchem::parser::pdb_read_options::all_models);
        lpdb.parse_molecule(ligand_mols);

        if (receptor_mols.size() == 1) {
            constant_receptor = true;
        } else if (ligand_mols.size() != receptor_mols.size()) {
            throw std::length_error("Differing number of complexes!");
        }
    }

    if (receptor_mols.get_idatm_types().size() == 1) {
        receptor_mols.compute_idatm_type()
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

    if (ligand_mols.get_idatm_types().size() == 1) {
        ligand_mols.compute_idatm_type()
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

    statchem::OMMIface::ForceField ffield;

    ffield.parse_gaff_dat_file(vm["gaff_dat"].as<std::string>())
        .parse_forcefield_file(vm["amber_xml"].as<std::string>())
        .parse_forcefield_file(vm["water_xml"].as<std::string>());

    auto mini_tol = vm["mini_tol"].as<double>();
    auto max_iter = vm["max_iter"].as<int>();

    auto platform = vm["platform"].as<std::string>();
    auto precision = vm["precision"].as<std::string>();
    auto accelerators = vm["accelerators"].as<std::string>();

    statchem::OMMIface::SystemTopology::loadPlugins();

    for (size_t i = 0; i < receptor_mols.size(); ++i) {
        statchem::molib::Molecule& protein =
            constant_receptor ? receptor_mols[0] : receptor_mols[i];
        statchem::molib::Molecule& ligand = ligand_mols[i];

        statchem::molib::Atom::Grid gridrec(protein.get_atoms());
        protein.prepare_for_mm(ffield, gridrec);

        ffield.insert_topology(protein);
        ffield.insert_topology(ligand);

        statchem::OMMIface::Modeler modeler(ffield, "phy", mini_tol, max_iter);

        modeler.add_topology(protein.get_atoms());
        modeler.add_topology(ligand.get_atoms());

        modeler.init_openmm(platform, precision, accelerators);

        modeler.add_crds(protein.get_atoms(), protein.get_crds());
        modeler.add_crds(ligand.get_atoms(), ligand.get_crds());

        modeler.unmask(ligand.get_atoms());
        modeler.unmask(protein.get_atoms());

        modeler.init_openmm_positions();

        modeler.minimize_state();

        // init with minimized coordinates
        statchem::molib::Molecule minimized_receptor(
            protein, modeler.get_state(protein.get_atoms()));
        statchem::molib::Molecule minimized_ligand(
            ligand, modeler.get_state(ligand.get_atoms()));

        minimized_receptor.undo_mm_specific();

        statchem::fileio::print_complex_pdb(std::cout, minimized_ligand,
                                            minimized_receptor, 0.000);

        ffield.erase_topology(ligand);
    }

    return 0;
}
