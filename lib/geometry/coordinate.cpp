/* This is coordinate.cpp and is part of StatChemLIB
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

#include "statchem/geometry/coordinate.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "statchem/geometry/matrix.hpp"
using namespace std;

namespace statchem {
namespace geometry {

unique_ptr<Coordinate> Coordinate::rotate(const Matrix& matrix) {
    unique_ptr<Coordinate> c(new Coordinate(*this));
    matrix.rotate(*c);
    return c;
}

unique_ptr<Coordinate> Coordinate::inverse_rotate(const Matrix& matrix) {
    unique_ptr<Coordinate> c(new Coordinate(*this));
    matrix.inverse_rotate(*c);
    return c;
}
}
}
