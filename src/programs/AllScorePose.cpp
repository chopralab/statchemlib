#include "AllScorePose.hpp"

#include <thread>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "statchem/fileio/inout.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
#include "statchem/score/score.hpp"

using namespace statchem_prog;

template <>
ProgramInfo statchem_prog::program_information<AllScorePose>() {
    return ProgramInfo("all_score_pose")
        .description(
            "Calculate the score of a given complex using all scoring "
            "functions.");
}

AllScorePose::AllScorePose() {}

int AllScorePose::run(int argc, char* argv[]) {
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

    po::options_description scoring_options("Scoring Function Arguments");
    scoring_options.add_options()(
        "dist", po::value<std::string>()->default_value(
                    "data/csd_complete_distance_distributions.txt.xz"),
        "Select one of the interatomic distance distribution file(s) "
        "provided with this program");

    po::options_description cmdln_options;
    cmdln_options.add(starting_inputs);
    cmdln_options.add(scoring_options);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdln_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout
            << "This program is used to calculate the KB-Score of a \n"
               "given compelx. To use it, specify a 'receptor' and a ligand'\n"
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

    auto dist = vm["dist"].as<std::string>();
    if (statchem::fileio::file_size(dist) <= 0) {
        throw std::runtime_error(std::string("File: '") + dist +
                                 std::string("' is invalid"));
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
            .erase_hydrogen();
    }

    if (ligand_mols.get_idatm_types().size() == 1) {
        ligand_mols.compute_idatm_type()
            .compute_hydrogen()
            .compute_bond_order()
            .compute_bond_gaff_type()
            .refine_idatm_type()
            .erase_hydrogen();
    }

    std::vector<std::string> distributions_file_raw;
    statchem::fileio::read_file(
        dist,
        distributions_file_raw);

    std::vector<std::string> scoring_names;
    std::map<std::string, std::unique_ptr<statchem::score::Score>> scoring_map;

    for (std::string comp : {"mean", "cumulative"}) {
        for (std::string func : {"radial", "normalized_frequency"}) {
            for (std::string ref : {"reduced", "complete"}) {
                for (auto cutoff : {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}) {
                    scoring_names.emplace_back(func + "_" + comp + "_" + ref +
                                               "_" + std::to_string(cutoff));
                    scoring_map[scoring_names.back()] =
                        std::unique_ptr<statchem::score::Score>(
                            new statchem::score::Score(comp, ref, func,
                                                       cutoff));

                    scoring_map[scoring_names.back()]
                        ->define_composition(receptor_mols.get_idatm_types(),
                                             ligand_mols.get_idatm_types())
                        .process_distributions(distributions_file_raw)
                        .compile_scoring_function();
                }
            }
        }
    }

    size_t num_threads = vm["ncpu"].as<int>() <= 0
                             ? std::thread::hardware_concurrency()
                             : static_cast<size_t>(vm["ncpu"].as<int>());

    std::vector<std::vector<double>> output(ligand_mols.size());

    std::vector<std::thread> threads;
    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        threads.push_back(std::thread([&, thread_id] {
            for (size_t i = thread_id; i < ligand_mols.size();
                 i += num_threads) {
                const auto& protein =
                    constant_receptor ? receptor_mols[0] : receptor_mols[i];
                const auto& ligand = ligand_mols[i];

                statchem::molib::Atom::Grid gridrec(protein.get_atoms());

                output[i].reserve(96);

                for (const auto& score_name : scoring_names) {
                    double new_score =
                        scoring_map[score_name]->non_bonded_energy(gridrec,
                                                                   ligand);
                    output[i].push_back(new_score);
                }
            }
        }));
    }

    for (auto&& thread : threads) {
        thread.join();
    }

    for (size_t i = 0; i < ligand_mols.size(); ++i) {
        std::cout << ligand_mols[i].name();
        for (auto score : output[i]) std::cout << ',' << score;
        std::cout << "\n";
    }

    return 0;
}
