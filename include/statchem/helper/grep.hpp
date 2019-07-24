/* This is grep.hpp and is part of StatChemLIB
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

#ifndef GREP_HPP
#define GREP_HPP

#include <istream>
#include <regex>

namespace statchem {

namespace grep {
// TODO: Make more like other inout functions???????
size_t count_matches(std::istream& stream, const std::regex& regex);
size_t find_first_case_greater_than(std::istream& stream,
                                    const std::regex& regex, double limit);
std::vector<std::string> search_stream(std::istream& stream,
                                       const std::regex& regex);
}
}

#endif
