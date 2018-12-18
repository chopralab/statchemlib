#include "PhysMinimize.hpp"

#include <iostream>
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

#include "programs/common.hpp"

using namespace statchem_prog;

template <>
ProgramInfo statchem_prog::program_information<PhysMinimize>() {
    return ProgramInfo("phys_min")
        .description("Minimize a complex using a physics-based potential.");
}

PhysMinimize::PhysMinimize() {}

bool PhysMinimize::process_options(int argc, char* argv[]) {
    auto starting_inputs = common_starting_inputs();

    auto ff_min = forcefield_options();
    auto openmm = openmm_options();

    po::options_description cmdln_options;
    cmdln_options.add(starting_inputs);
    cmdln_options.add(ff_min);
    cmdln_options.add(openmm);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdln_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        __help_text << "This program performs energy minization using"
                    << "the Amber forcefield.\n\n";
        __help_text << starting_inputs_help() << "\n";
        __help_text << cmdln_options << std::endl;
        return false;
    }

    __constant_receptor =
        process_starting_inputs(vm, __receptor_mols, __ligand_mols);

    process_forcefield_options(vm, __ffield, __mini_tol, __iter_max);
    process_openmm_options(vm, __platform, __precision, __accelerators, __checkpoint);

    return true;
}

int PhysMinimize::run() {
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

    statchem::OMMIface::SystemTopology::loadPlugins();

    for (size_t i = 0; i < __ligand_mols.size(); ++i) {
        statchem::molib::Molecule& protein =
            __constant_receptor ? __receptor_mols[0] : __receptor_mols[i];
        statchem::molib::Molecule& ligand = __ligand_mols[i];

        statchem::molib::Atom::Grid gridrec(protein.get_atoms());
        protein.prepare_for_mm(__ffield, gridrec);

        __ffield.insert_topology(protein);
        __ffield.insert_topology(ligand);

        statchem::OMMIface::Modeler modeler(__ffield, "phy", 0.0,
                                            __mini_tol, __iter_max);

        modeler.add_topology(protein.get_atoms());
        modeler.add_topology(ligand.get_atoms());

        modeler.init_openmm(__platform, __precision, __accelerators);

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

        __ffield.erase_topology(ligand);
    }

    return 0;
}
