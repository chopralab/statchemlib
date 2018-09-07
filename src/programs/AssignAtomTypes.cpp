#include "AssignAtomTypes.hpp"

#include <iostream>
#include <thread>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "statchem/fileio/inout.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
#include "statchem/score/score.hpp"

#include "programs/common.hpp"

using namespace statchem_prog;

template <>
ProgramInfo statchem_prog::program_information<AssignAtomTypes>() {
    return ProgramInfo("assign_atom_types")
        .description(
            "Determine atom types for a molecule.");
}

AssignAtomTypes::AssignAtomTypes() {}

bool AssignAtomTypes::process_options(int argc, char* argv[]) {
    po::options_description starting_inputs("File input/output options");
    starting_inputs.add_options()("help,h", "Show this help menu.")(
        "ligand,l", po::value<std::string>()->default_value("ligand.mol2"),
        "Ligand filename. This can be in either PDB or MOL2 format.")(
        "output,o", po::value<std::string>()->default_value("output.pdb"),
        "Output filename. Must be in the PDB format due to custom types");

    po::options_description cmdln_options;
    cmdln_options.add(starting_inputs);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdln_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        __help_text << "This program assigns atom types to a molecule "
                    << "for use in functions provided by StatChemLib.\n\n";
        __help_text << cmdln_options << std::endl;
        return false;
    }

    auto ligand = vm["ligand"].as<std::string>();

    statchem::parser::FileParser lpdb(
        ligand, statchem::parser::pdb_read_options::all_models | statchem::parser::pdb_read_options::hydrogens);
    lpdb.parse_molecule(__molecules);

    __output = vm["output"].as<std::string>();

    return true;
}


int AssignAtomTypes::run() {
    __molecules.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .compute_hydrogen()
        .compute_ring_type()
        .compute_gaff_type()
        .erase_temporary_hydrogen();

    std::cout << __molecules;

    return 0;
}
