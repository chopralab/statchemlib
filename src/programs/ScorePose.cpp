#include "ScorePose.hpp"

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
ProgramInfo statchem_prog::program_information<ScorePose>() {
    return ProgramInfo("score_pose")
        .description(
            "Calculate the score of a given complex using a single scoring "
            "function.");
}

ScorePose::ScorePose() {}

bool ScorePose::process_options(int argc, char* argv[]) {
    auto starting_inputs = common_starting_inputs();

    auto score_options = scoring_options();

    po::options_description cmdln_options;
    cmdln_options.add(starting_inputs);
    cmdln_options.add(score_options);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdln_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        __help_text << "This program computes the score of a complex using "
                    << " a specified scoring function.\n\n";
        __help_text << starting_inputs_help() << "\n";
        __help_text << cmdln_options << std::endl;
        return false;
    }

    __constant_receptor =
        process_starting_inputs(vm, __receptor_mols, __ligand_mols);

    __num_threads = vm["ncpu"].as<int>() <= 0
                        ? std::thread::hardware_concurrency()
                        : static_cast<size_t>(vm["ncpu"].as<int>());

    process_scorint_options(vm, __ref, __comp, __func, __cutoff, __dist);

    return true;
}

int ScorePose::run() {
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

    __score = std::unique_ptr<statchem::score::Score>(
        new statchem::score::Score(__ref, __comp, __func, __cutoff));

    __score
        ->define_composition(__receptor_mols.get_idatm_types(),
                             __ligand_mols.get_idatm_types())
        .process_distributions(__dist)
        .compile_scoring_function();

    std::vector<double> output(__ligand_mols.size());

    std::vector<std::thread> threads;
    for (size_t thread_id = 0; thread_id < __num_threads; ++thread_id) {
        threads.push_back(std::thread([&, thread_id] {
            for (size_t i = thread_id; i < __ligand_mols.size();
                 i += __num_threads) {
                const auto& protein = __constant_receptor ? __receptor_mols[0]
                                                          : __receptor_mols[i];
                const auto& ligand = __ligand_mols[i];

                statchem::molib::Atom::Grid gridrec(protein.get_atoms());

                output[i] = __score->non_bonded_energy(gridrec, ligand);
            }
        }));
    }

    for (auto&& thread : threads) {
        thread.join();
    }

    for (size_t i = 0; i < output.size(); ++i) {
        std::cout << output[i] << '\n';
    }

    return 0;
}
