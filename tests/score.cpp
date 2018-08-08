#include "statchem/parser/fileparser.hpp"
#include "statchem/helper/grep.hpp"
#include "statchem/score/score.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("RMR6 Score") {
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
        .erase_hydrogen()    // needed because refine changes connectivities
        .compute_hydrogen()  // needed because refine changes connectivities
        .compute_ring_type()
        .compute_gaff_type()
        .compute_rotatable_bonds()
        .erase_hydrogen();

    rmol.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen()    // needed because refine changes connectivities
        .compute_hydrogen()  // needed because refine changes connectivities
        .compute_ring_type()
        .compute_gaff_type()
        .compute_rotatable_bonds()
        .erase_hydrogen();

    statchem::score::Score score( "mean", "reduced", "radial", 6);
    score
        .define_composition(rmol.get_idatm_types(),
                            lmol.get_idatm_types())
        .process_distributions_file("../data/csd_complete_distance_distributions.txt.xz")
        .compile_scoring_function();

    statchem::molib::Atom::Grid gridrec(rmol[0].get_atoms());
    auto my_score = score.non_bonded_energy(gridrec, lmol[0]);
    CHECK(std::fabs(my_score - (-4.832775)) < 1e-6 );
}
