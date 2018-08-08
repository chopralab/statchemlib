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
