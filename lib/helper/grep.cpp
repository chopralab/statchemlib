/* This is grep.cpp and is part of StatChemLIB
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

#include "statchem/helper/grep.hpp"

namespace statchem {
size_t grep::count_matches(std::istream& stream, const std::regex& regex) {
    size_t count = 0;

    std::smatch what;
    std::string line;

    while (std::getline(stream, line)) {
        bool result = std::regex_search(line, what, regex);
        if (result) {
            ++count;
        }
    }

    return count;
}

size_t grep::find_first_case_greater_than(std::istream& stream,
                                          const std::regex& regex,
                                          double limit) {
    size_t count = 0;

    std::smatch what;
    std::string line;

    while (std::getline(stream, line)) {
        bool result = std::regex_match(line, what, regex);
        if (result) {
            std::string result(what[1].first, what[1].second);
            double value = std::stod(result);

            if (value > limit) {
                break;
            }

            ++count;
        }
    }

    return count;
}

std::vector<std::string> grep::search_stream(std::istream& stream,
                                             const std::regex& regex) {
    std::vector<std::string> all_matches;

    std::match_results<std::string::const_iterator> what;

    std::string line;

    while (std::getline(stream, line)) {
        bool result = std::regex_search(line, what, regex);
        if (result) {
            all_matches.push_back(std::string(what[1].first, what[1].second));
        }
    }

    return all_matches;
}
}
