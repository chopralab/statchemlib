#include "statchem/parser/fileparser.hpp"
#include "statchem/helper/grep.hpp"
#include "statchem/helper/help.hpp"

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Read a mol2 file and assign rotable bonds") {
    statchem::parser::FileParser lmol2("files/drugs.mol2");
    statchem::molib::Molecules mols;

    lmol2.parse_molecule(mols);
    CHECK(mols.size() == 3);
    CHECK(mols[0].name() == "random");
    CHECK(mols[1].name() == "tibolone");
    CHECK(mols[2].name() == "loratadine");

    CHECK(mols[0].first().first().first().first().size() == 28);
    CHECK(mols[1].first().first().first().first().size() == 23);

    mols.compute_idatm_type()
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

    CHECK(mols[0].first().first().first().first().size() == 28);
    CHECK(mols[1].first().first().first().first().size() == 23);

    auto atomtypes = mols.get_idatm_types();
    CHECK(atomtypes.size() == 14);
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("C1")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("C2")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("C3")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("Car")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("N2")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("N3+")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("Npl")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("O2")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("O3")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("O3-")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("Oar")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("F")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("Cl")));
    CHECK(atomtypes.count(statchem::help::idatm_mask.at("Pac")));

    auto path = fs::temp_directory_path() / fs::unique_path();

    mols.compute_overlapping_rigid_segments(path.string());
    std::set<int> added;
    statchem::molib::Molecules frags;
    statchem::molib::create_mols_from_seeds(added, frags, mols);

    CHECK(frags.size() == added.size());

    std::regex reg("$");
    std::ifstream seeds_file(path.string());
    auto count = statchem::grep::count_matches(seeds_file, reg);

    CHECK(count == frags.size());
    fs::remove(path);
}
