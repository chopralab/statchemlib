#ifndef _STCH_COMMON_OPTIONS_HPP_
#define _STCH_COMMON_OPTIONS_HPP_

#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include "statchem/fileio/inout.hpp"
#include "statchem/modeler/forcefield.hpp"
#include "statchem/parser/fileparser.hpp"

namespace statchem_prog {
namespace po = boost::program_options;

inline po::options_description common_starting_inputs() {
    po::options_description starting_inputs("Starting input files");
    starting_inputs.add_options()("help,h", "Show this help menu.")(
        "receptor,r", po::value<std::string>()->default_value("receptor.pdb"),
        "Receptor filename. Will be ignored if the complex option is given")(
        "ligand,l", po::value<std::string>()->default_value("ligand.mol2"),
        "Ligand filename. Will be ignored if the complex option is given")(
        "complex,c", po::value<std::string>(),
        "File with both the ligand and receptor present. Overrides above "
        "options")(
        "ncpu,n", po::value<int>()->default_value(-1),
        "Number of CPUs to use concurrently (use -1 to use all CPUs)");

    return starting_inputs;
}

inline std::string starting_inputs_help() {
    return "This program operates on one or many combinations of a 'receptor'\n"
           "and a 'ligand'. There are few methods for specifying a receptor\n"
           "and a ligand' file. First, you can specify a complex file in the\n"
           "PDB format which separates ligand from repector with TER records.\n"
           "You may also specify two files containing the components. The\n"
           "following cases are valid:\n\n"
           "(1) A single complex PDB file with multiple complexes "
           "separted by MODEL and ENDMDL.\n"
           "\tThe ligand is separated with a TER record.\n"
           "(2) A receptor file containing a single biomolecule and "
           "a ligand file containing\n\tmutliple small-molecules\n"
           "(3) A receptor file containing multiple biomolecules and "
           "a ligand file containing\n\tthe same number of "
           "small-molecules\n";
}

inline bool process_starting_inputs(po::variables_map& vm,
                                    statchem::molib::Molecules& rec_mols,
                                    statchem::molib::Molecules& lig_mols) {
    if (vm.count("complex")) {
        auto complex = vm["complex"].as<std::string>();

        statchem::parser::FileParser drpdb(
            complex, statchem::parser::pdb_read_options::protein_poses_only |
                         statchem::parser::pdb_read_options::all_models);

        drpdb.parse_molecule(rec_mols);

        statchem::parser::FileParser dlpdb(
            complex, statchem::parser::pdb_read_options::docked_poses_only |
                         statchem::parser::pdb_read_options::skip_atom |
                         statchem::parser::pdb_read_options::all_models);
        dlpdb.parse_molecule(lig_mols);

    } else {
        auto receptor = vm["receptor"].as<std::string>();
        statchem::parser::FileParser rpdb(
            receptor, statchem::parser::pdb_read_options::all_models);
        rpdb.parse_molecule(rec_mols);

        auto ligand = vm["ligand"].as<std::string>();
        statchem::parser::FileParser lpdb(
            ligand, statchem::parser::pdb_read_options::all_models);
        lpdb.parse_molecule(lig_mols);

        if (rec_mols.size() == 1) {
            return true;
        } else if (lig_mols.size() != rec_mols.size()) {
            throw std::length_error("Differing number of complexes!");
        }
    }

    return false;
}

inline po::options_description forcefield_options() {
    po::options_description ff_min("Forcefield and Minimization Options");
    ff_min.add_options()(
        "amber_xml",
        po::value<std::string>()->default_value("data/amber10.xml"),
        "Receptor XML parameters (and topology) input file")(
        "water_xml", po::value<std::string>()->default_value("data/tip3p.xml"),
        "Water XML parameters (and topology) input file")(
        "gaff_dat", po::value<std::string>()->default_value("data/gaff.dat"),
        "Gaff DAT forcefield input file")(
        "mini_tol", po::value<double>()->default_value(0.0001),
        "Minimization tolerance")(
        "max_iter", po::value<int>()->default_value(10),
        "Maximum iterations for minimization during linking");

    return ff_min;
}

inline void process_forcefield_options(po::variables_map& vm,
                                       statchem::OMMIface::ForceField& ffield,
                                       double& mini_tol, int& max_iter) {
    ffield.parse_gaff_dat_file(vm["gaff_dat"].as<std::string>())
        .parse_forcefield_file(vm["amber_xml"].as<std::string>())
        .parse_forcefield_file(vm["water_xml"].as<std::string>());

    mini_tol = vm["mini_tol"].as<double>();
    if (mini_tol < 0.0) {
        throw po::validation_error(po::validation_error::invalid_option_value,
                                   "mini_tol", std::to_string(mini_tol));
    }

    max_iter = vm["max_iter"].as<int>();
    if (max_iter < 0.0) {
        throw po::validation_error(po::validation_error::invalid_option_value,
                                   "max_iter", std::to_string(max_iter));
    }
}

inline po::options_description openmm_options() {
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

    return openmm;
}

inline void process_openmm_options(po::variables_map& vm, std::string& platform,
                                   std::string& precision,
                                   std::string& accelerators) {
    platform = vm["platform"].as<std::string>();
    precision = vm["precision"].as<std::string>();
    accelerators = vm["accelerators"].as<std::string>();
}

inline po::options_description scoring_options() {
    po::options_description scoring_options("Scoring Function Arguments");
    scoring_options.add_options()(
        "dist",
        po::value<std::string>()->default_value(
            "data/csd_complete_distance_distributions.txt.xz"),
        "Select one of the interatomic distance distribution file(s) "
        "provided with this program")(
        "ref", po::value<std::string>()->default_value("mean"),
        "Normalization method for the reference state ('mean' is averaged "
        "over all atom type pairs, whereas 'cumulative' is a summation for "
        "atom type pairs)")(
        "comp", po::value<std::string>()->default_value("reduced"),
        "Atom types used in calculating reference state 'reduced' or "
        "'complete'"
        "('reduced' includes only those atom types present in the "
        "specified receptor and small molecule, whereas 'complete' "
        "includes all atom types)")(
        "func", po::value<std::string>()->default_value("radial"),
        "Function for calculating scores 'radial' or "
        "'normalized_frequency'")("cutoff", po::value<int>()->default_value(6),
                                  "Cutoff length [4-15].");

    return scoring_options;
}

inline void process_scorint_options(po::variables_map& vm, std::string& ref,
                                    std::string& comp, std::string& func,
                                    double& cutoff, std::string& dist) {
    ref = vm["ref"].as<std::string>();
    comp = vm["comp"].as<std::string>();
    func = vm["func"].as<std::string>();
    int cutoff_temp = vm["cutoff"].as<int>();

    if (ref != "mean" && ref != "cumulative") {
        throw po::validation_error(po::validation_error::invalid_option_value,
                                   "ref", ref);
    }

    if (comp != "reduced" && comp != "complete") {
        throw po::validation_error(po::validation_error::invalid_option_value,
                                   "comp", comp);
    }

    if (func != "radial" && func != "normalized_frequency") {
        throw po::validation_error(po::validation_error::invalid_option_value,
                                   "func", func);
    }

    if (cutoff_temp < 4 || cutoff_temp > 15) {
        throw po::validation_error(po::validation_error::invalid_option_value,
                                   "cutoff", std::to_string(cutoff));
    }

    cutoff = static_cast<double>(cutoff_temp);

    dist = vm["dist"].as<std::string>();
    if (statchem::fileio::file_size(dist) <= 0) {
        throw std::runtime_error(std::string("File: '") + dist +
                                 std::string("' is invalid"));
    }
}

}  // namespace statchem_prog

#endif
