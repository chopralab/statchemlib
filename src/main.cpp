#include <iomanip>
#include <iostream>
#include <string>

#include "ProgramManager.hpp"

using namespace statchem_prog;

int main(int argc, char* argv[]) {
    if (argc <= 1) {
        std::cerr << "You must specify a program to run." << std::endl;
        std::cerr << "Run '" << argv[0] << " help' to see availible programs."
                  << std::endl;
        return 1;
    }

    std::string program_name(argv[1]);
    if (program_name == "help" || program_name == "-h" ||
        program_name == "--help") {
        std::cout << "Avaible programs are:\n" << std::endl;
        auto programs = ProgramManager::get().programs();
        for (auto prog : programs) {
            std::cout << "    " << std::left << std::setw(15) << prog.name();
            std::cout << "\t" << prog.description() << std::endl;
        }
        std::cout << std::endl
                  << "To see help for a specific program, run" << std::endl;
        std::cout << argv[0] << " <program_name> --help" << std::endl;
        return 0;
    }

    try {
        auto program = ProgramManager::get().name(argv[1]);

        return program()->run(argc - 1, argv + 1);
    } catch(std::exception& e) {
        std::cerr << "'" << argv[0] << " " << argv[1] << "' failed because: " << std::endl;
        std::cerr << e.what() << std::endl;
        return 2;
    }
}
