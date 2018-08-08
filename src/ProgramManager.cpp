// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cctype>
#include <sstream>

#include "ProgramManager.hpp"

#include "programs/ScorePose.hpp"

using namespace statchem_prog;

ProgramManager::ProgramManager() {
    // Add formats here
    this->add_format<ScorePose>();
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
