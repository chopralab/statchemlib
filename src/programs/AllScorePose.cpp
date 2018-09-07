#include "AllScorePose.hpp"

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
ProgramInfo statchem_prog::program_information<AllScorePose>() {
    return ProgramInfo("all_score_pose")
        .description(
            "Calculate the score of a given complex using all scoring "
            "functions.");
}

AllScorePose::AllScorePose() {}

bool AllScorePose::process_options(int argc, char* argv[]) {
    auto starting_inputs = common_starting_inputs();
    starting_inputs.add_options()(
        "ncpu,n", po::value<int>()->default_value(-1),
        "Number of CPUs to use concurrently (use -1 to use all CPUs)");

    po::options_description scoring_options("Scoring Function Arguments");
    scoring_options.add_options()(
        "dist",
        po::value<std::string>(&__dist)->default_value(
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
        __help_text << "This program scores complexes using all the scoring "
                    << "functions provided by StatChemLib.\n\n";
        __help_text << starting_inputs_help() << "\n";
        __help_text << cmdln_options << std::endl;
        return false;
    }

    __constant_receptor =
        process_starting_inputs(vm, __receptor_mols, __ligand_mols);

    __num_threads = vm["ncpu"].as<int>() <= 0
                        ? std::thread::hardware_concurrency()
                        : static_cast<size_t>(vm["ncpu"].as<int>());

    return true;
}

int AllScorePose::run() {
    if (__receptor_mols.get_idatm_types().size() == 1) {
        __receptor_mols.compute_idatm_type()
            .compute_hydrogen()
            .compute_bond_order()
            .compute_bond_gaff_type()
            .refine_idatm_type()
            .erase_hydrogen();
    }

    if (__ligand_mols.get_idatm_types().size() == 1) {
        __ligand_mols.compute_idatm_type()
            .compute_hydrogen()
            .compute_bond_order()
            .compute_bond_gaff_type()
            .refine_idatm_type()
            .erase_hydrogen();
    }

    statchem::score::AtomicDistributions distributions(__dist);

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
                        ->define_composition(__receptor_mols.get_idatm_types(),
                                             __ligand_mols.get_idatm_types())
                        .process_distributions(distributions)
                        .compile_scoring_function();
                }
            }
        }
    }

    std::vector<std::vector<double>> output(__ligand_mols.size());

    std::vector<std::thread> threads;
    for (size_t thread_id = 0; thread_id < __num_threads; ++thread_id) {
        threads.push_back(std::thread([&, thread_id] {
            for (size_t i = thread_id; i < __ligand_mols.size();
                 i += __num_threads) {
                const auto& protein = __constant_receptor ? __receptor_mols[0]
                                                          : __receptor_mols[i];
                const auto& ligand = __ligand_mols[i];

                const statchem::molib::Atom::Grid gridrec(protein.get_atoms());

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

    for (size_t i = 0; i < __ligand_mols.size(); ++i) {
        std::cout << __ligand_mols[i].name();
        for (auto score : output[i]) std::cout << ',' << score;
        std::cout << "\n";
    }

    return 0;
}
