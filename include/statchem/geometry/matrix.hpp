/* This is matrix.hpp and is part of StatChemLIB
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

#ifndef MATRIX_H
#define MATRIX_H
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include "statchem/helper/error.hpp"

namespace statchem {

namespace geometry {

class Coordinate;

class Matrix {
    struct MatrixPrivate;
    std::unique_ptr<MatrixPrivate> __matrix;

   public:
    typedef std::tuple<double, double, double, double> matrix_tuple;
    Matrix();
    Matrix(const double d[9]);
    Matrix(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    ~Matrix();
    void set_row(const int& row, const matrix_tuple& t);

    void rotate(Coordinate& coord) const;
    void inverse_rotate(Coordinate& coord) const;

    friend std::ostream& operator<<(std::ostream& stream, const Matrix& m);
};
}
}

#endif
