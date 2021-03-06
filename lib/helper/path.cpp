/* This is path.cpp and is part of StatChemLIB
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

#include "statchem/helper/path.hpp"
#include <math.h>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>

namespace statchem {
std::string Path::join(const std::string& str1, const std::string& str2) {
    boost::filesystem::path p1(str1);
    boost::filesystem::path p2(str2);
    p1 /= p2;
    return p1.string();
}
}
