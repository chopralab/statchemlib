#include "MakeObjective.hpp"

#include <thread>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "statchem/fileio/inout.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
#include "statchem/score/score.hpp"

#include "programs/common.hpp"

using namespace statchem_prog;

template <>
ProgramInfo statchem_prog::program_information<MakeObjective>() {
    return ProgramInfo("make_objective")
        .description(
            "Generate an objective function based on a statistical potential.");
}

MakeObjective::MakeObjective() {}

bool MakeObjective::process_options(int argc, char* argv[]) {
    auto starting_inputs = common_starting_inputs();

    auto score_options = scoring_options();
    score_options.add_options()(
        "step_size", po::value<double>(&__step_size)->default_value(0.01),
        "Step size for fitted bsplines.");

    std::vector<std::string> atom_type_names;
    auto objective_options = po::options_description("Objective Options");
    objective_options.add_options()(
        "atom_types",
        po::value<std::vector<std::string>>(&atom_type_names)->multitoken(),
        "Which atom types to use for compiling the objective function.")(
        "obj_dir", po::value<std::string>(&__obj_dir)->default_value("obj"),
        "Location to store output files");

    po::options_description cmdln_options;
    cmdln_options.add(starting_inputs);
    cmdln_options.add(score_options);
    cmdln_options.add(objective_options);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdln_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        __help_text << "This program compiles an objective function using "
                    << " a specified scoring function.\n\n";
        __help_text << cmdln_options << std::endl;
        return false;
    }

    process_scoring_options(vm, __ref, __comp, __func, __cutoff, __dist);

    if (!atom_type_names.empty()) {
        for (auto a : atom_type_names) {
            auto pos = statchem::help::idatm_mask.find(a);

            if (pos == statchem::help::idatm_mask.end()) {
                throw po::validation_error(
                    po::validation_error::invalid_option_value, "atom_types",
                    a);
            }

            __atom_types.insert(pos->second);
        }
    } else {
        for (auto b : statchem::help::idatm_mask) {
            __atom_types.insert(b.second);
        }
    }

    return true;
}

int MakeObjective::run() {
    __score = std::unique_ptr<statchem::score::KBFF>(new statchem::score::KBFF(
        __ref, __comp, __func, __cutoff, __step_size));

    __score->define_composition(__atom_types, __atom_types)
        .process_distributions(__dist)
        .compile_scoring_function();
    __score->compile_objective_function();

    __score->output_objective_function(__obj_dir);

    return 0;
}
