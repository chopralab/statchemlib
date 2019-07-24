/* This is ProgramManager.cpp and is part of StatChemLIB
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

#include <cctype>
#include <sstream>

#include "ProgramManager.hpp"

#include "programs/ScorePose.hpp"
#include "programs/AllScorePose.hpp"
#include "programs/KBMinimize.hpp"
#include "programs/PhysMinimize.hpp"
#include "programs/MakeObjective.hpp"
#include "programs/KBDynamics.hpp"
#include "programs/PhysDynamics.hpp"
#include "programs/AssignAtomTypes.hpp"

using namespace statchem_prog;

ProgramManager::ProgramManager() {
    // Add formats here
    this->add_format<ScorePose>();
    this->add_format<AllScorePose>();
    this->add_format<KBMinimize>();
    this->add_format<PhysMinimize>();
    this->add_format<MakeObjective>();
    this->add_format<KBDynamics>();
    this->add_format<PhysDynamics>();
    this->add_format<AssignAtomTypes>();
}

ProgramManager& ProgramManager::get() {
    static ProgramManager _instance;
    return _instance;
}

ProgramManager::iterator ProgramManager::find_name(const program_map_t& formats,
                                                   const std::string& name) {
    for (auto it = formats.begin(); it != formats.end(); it++) {
        if (it->first.name() == name) {
            return it;
        }
    }
    return formats.end();
}

void ProgramManager::register_program(ProgramInfo info,
                                     program_creator_t creator) {
    if (info.name() == "") {
        throw std::logic_error("can not program a format with no name");
    } else if (info.name() != "" &&
               find_name(__programs, info.name()) != __programs.end()) {
        throw std::domain_error("Attempt to add a format twice!");
    } else {
        __programs.push_back({info, creator});
    }
}

program_creator_t ProgramManager::name(const std::string& name) {
    auto it = find_name(__programs, name);
    if (it == __programs.end()) {
        throw std::runtime_error(
            "Cannot find the specificed program. Check the 'help' to see "
            "availible options");
    }
    return it->second;
}

std::vector<ProgramInfo> ProgramManager::programs() {
    auto metadata = std::vector<ProgramInfo>();
    metadata.reserve(__programs.size());
    for (auto& pair : __programs) {
        metadata.emplace_back(pair.first);
    }
    return metadata;
}
