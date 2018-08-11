#include "statchem/score/score.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/grep.hpp"
#include "statchem/parser/fileparser.hpp"

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
        .erase_hydrogen();

    rmol.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen();

    statchem::score::Score score("mean", "reduced", "radial", 6);
    score.define_composition(rmol.get_idatm_types(), lmol.get_idatm_types())
        .process_distributions(
            "../data/csd_complete_distance_distributions.txt.xz")
        .compile_scoring_function();

    statchem::molib::Atom::Grid gridrec(rmol[0].get_atoms());
    auto my_score = score.non_bonded_energy(gridrec, lmol[0]);
    std::cout << my_score << std::endl;
    CHECK(std::fabs(my_score - (-4.832775)) < 1e-6);
}

TEST_CASE("Muliple Scoring Functions") {
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
        .erase_hydrogen();

    rmol.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen();

    statchem::score::AtomicDistributions distributions(
        "../data/csd_complete_distance_distributions.txt.xz");
    auto ltypes = lmol.get_idatm_types();
    auto rtypes = rmol.get_idatm_types();

    statchem::score::Score rmc15("mean", "complete", "radial", 15);
    rmc15.define_composition(rtypes, ltypes)
        .process_distributions(distributions)
        .compile_scoring_function();

    CHECK_THROWS(rmc15.define_composition(rtypes, ltypes));
    CHECK_THROWS(rmc15.process_distributions(distributions));

    statchem::score::Score fmc10("mean", "complete", "normalized_frequency",
                                 10);
    fmc10.define_composition(rtypes, ltypes)
        .process_distributions(distributions)
        .compile_scoring_function();

    statchem::score::Score fcc4("cumulative", "complete",
                                "normalized_frequency", 4);
    fcc4.define_composition(rtypes, ltypes)
        .process_distributions(distributions)
        .compile_scoring_function();

    statchem::molib::Atom::Grid gridrec(rmol[0].get_atoms());
    auto rmc15_score = rmc15.non_bonded_energy(gridrec, lmol[0]);
    auto fmc10_score = fmc10.non_bonded_energy(gridrec, lmol[0]);
    auto fcc4_score = fcc4.non_bonded_energy(gridrec, lmol[0]);
    CHECK(std::fabs(rmc15_score - (-10932.8691688)) < 1e-6);
    CHECK(std::fabs(fmc10_score - (-3962.8519988)) < 1e-6);
    CHECK(std::fabs(fcc4_score - (0.0451727)) < 1e-6);
}
