#include "statchem/modeler/forcefield.hpp"
#include "statchem/modeler/modeler.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
#include "statchem/score/kbff.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Knowledge-based energy minization") {
    statchem::parser::FileParser lpdb("files/6drw_lig.pdb");
    statchem::parser::FileParser rpdb("files/6drw.pdb");

    statchem::molib::Molecules lmol;
    statchem::molib::Molecules rmol;
    lpdb.parse_molecule(lmol);
    rpdb.parse_molecule(rmol);

    lmol.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen()
        .compute_hydrogen()
        .compute_ring_type()
        .compute_gaff_type()
        .erase_hydrogen();

    rmol.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen()
        .compute_hydrogen()
        .compute_ring_type()
        .compute_gaff_type()
        .erase_hydrogen();

    statchem::score::AtomicDistributions distributions(
        "../data/csd_complete_distance_distributions.txt.xz");

    statchem::score::Score scoring_func("mean", "reduced", "radial", 6);
    scoring_func
        .define_composition(rmol.get_idatm_types(), lmol.get_idatm_types())
        .process_distributions(distributions)
        .compile_scoring_function();

    statchem::score::KBFF objective_func("mean", "complete", "radial", 15,
                                         0.01);
    objective_func
        .define_composition(rmol.get_idatm_types(), lmol.get_idatm_types())
        .process_distributions(distributions)
        .compile_scoring_function();
    objective_func.compile_objective_function();

    statchem::OMMIface::ForceField ffield;

    ffield.parse_gaff_dat_file("../data/gaff.dat")
        .add_kb_forcefield(objective_func, 6.0)
        .parse_forcefield_file("../data/amber10.xml")
        .parse_forcefield_file("../data/tip3p.xml");

    statchem::molib::Molecule& protein = rmol[0];
    statchem::molib::Molecule& ligand = lmol[0];

    statchem::molib::Atom::Grid gridrec(protein.get_atoms());
    protein.prepare_for_mm(ffield, gridrec);

    ffield.insert_topology(protein);
    ffield.insert_topology(ligand);

    statchem::OMMIface::Modeler modeler(ffield, "kb", 1.0, 0.00001, 100);

    modeler.add_topology(protein.get_atoms());
    modeler.add_topology(ligand.get_atoms());

    modeler.init_openmm("Reference");

    modeler.add_crds(protein.get_atoms(), protein.get_crds());
    modeler.add_crds(ligand.get_atoms(), ligand.get_crds());

    modeler.unmask(ligand.get_atoms());
    modeler.unmask(protein.get_atoms());

    modeler.init_openmm_positions();
    double potential = modeler.potential_energy();
    modeler.minimize_state();

    // init with minimized coordinates
    statchem::molib::Molecule minimized_receptor(
        protein, modeler.get_state(protein.get_atoms()));
    statchem::molib::Molecule minimized_ligand(
        ligand, modeler.get_state(ligand.get_atoms()));

    minimized_receptor.undo_mm_specific();

    statchem::molib::Atom::Grid gridrec_min(minimized_receptor.get_atoms());
    double min_energy =
        scoring_func.non_bonded_energy(gridrec_min, minimized_ligand);
    double min_potential = modeler.potential_energy();

    CHECK(min_potential < potential);
    CHECK(std::fabs(min_energy - -5.51348) < 1e-3);
}

TEST_CASE("Physics-based energy minization") {
    statchem::parser::FileParser lpdb("files/6drw_lig.pdb");
    statchem::parser::FileParser rpdb("files/6drw.pdb");

    statchem::molib::Molecules lmol;
    statchem::molib::Molecules rmol;
    lpdb.parse_molecule(lmol);
    rpdb.parse_molecule(rmol);

    lmol.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen()
        .compute_hydrogen()
        .compute_ring_type()
        .compute_gaff_type()
        .erase_hydrogen();

    rmol.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen()
        .compute_hydrogen()
        .compute_ring_type()
        .compute_gaff_type()
        .erase_hydrogen();

    statchem::OMMIface::ForceField ffield;

    ffield.parse_gaff_dat_file("../data/gaff.dat")
        .parse_forcefield_file("../data/amber10.xml")
        .parse_forcefield_file("../data/tip3p.xml");

    statchem::molib::Molecule& protein = rmol[0];
    statchem::molib::Molecule& ligand = lmol[0];

    statchem::molib::Atom::Grid gridrec(protein.get_atoms());
    protein.prepare_for_mm(ffield, gridrec);

    ffield.insert_topology(protein);
    ffield.insert_topology(ligand);

    statchem::OMMIface::Modeler modeler(ffield, "phy", 1.0, 0.00001, 100);

    modeler.add_topology(protein.get_atoms());
    modeler.add_topology(ligand.get_atoms());

    modeler.init_openmm("Reference");

    modeler.add_crds(protein.get_atoms(), protein.get_crds());
    modeler.add_crds(ligand.get_atoms(), ligand.get_crds());

    modeler.unmask(ligand.get_atoms());
    modeler.unmask(protein.get_atoms());

    modeler.init_openmm_positions();
    double potential = modeler.potential_energy();
    modeler.minimize_state();

    // init with minimized coordinates
    statchem::molib::Molecule minimized_receptor(
        protein, modeler.get_state(protein.get_atoms()));
    statchem::molib::Molecule minimized_ligand(
        ligand, modeler.get_state(ligand.get_atoms()));

    minimized_receptor.undo_mm_specific();

    double min_potential = modeler.potential_energy();

    CHECK(min_potential < potential);
}
