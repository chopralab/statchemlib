/* This is main.cpp and is part of StatChemLIB
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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
        return 1;
    }

    try {
        auto program_factory = ProgramManager::get().name(program_name);
        auto program = program_factory();
        if (program->process_options(argc - 1, argv + 1)) {
            return program->run();
        }

        std::cout << program->get_help() << std::endl;
        return 1;
    } catch(std::exception& e) {
        std::cerr << "'" << argv[0] << " " << argv[1] << "' failed because: " << std::endl;
        std::cerr << e.what() << std::endl;
        std::cerr << "Please run '" << argv[0] << " " << argv[1] << " --help' for more information." << std::endl;
        return 2;
    }
}
