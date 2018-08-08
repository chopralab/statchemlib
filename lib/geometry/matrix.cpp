#include "statchem/geometry/matrix.hpp"
#include "statchem/geometry/coordinate.hpp"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector_double.h>

namespace statchem {
namespace geometry {

struct Matrix::MatrixPrivate {
    gsl_matrix* first;
    gsl_vector* second;
   public:
    MatrixPrivate() : first(gsl_matrix_alloc(3, 3)), second(gsl_vector_alloc(3)) {}
    ~MatrixPrivate() {
        gsl_matrix_free(first);
        gsl_vector_free(second);
    }
};

Matrix::Matrix() : __matrix(new MatrixPrivate) {}

Matrix::Matrix(const double d[9]) : __matrix(new MatrixPrivate) {
    for (int i = 0; i < 9; i++) {
        gsl_matrix_set(__matrix->first, i / 3, i % 3, d[i]);
    }
    gsl_vector_set_all(__matrix->second, 0);
}

Matrix::Matrix(const Matrix& other) : __matrix(new MatrixPrivate) {
    gsl_matrix_memcpy(__matrix->first, other.__matrix->first);
    gsl_vector_memcpy(__matrix->second, other.__matrix->second);
}

Matrix& Matrix::operator=(const Matrix& other) {
    gsl_matrix_memcpy(__matrix->first, other.__matrix->first);
    gsl_vector_memcpy(__matrix->second, other.__matrix->second);
    return *this;
}

Matrix::~Matrix() {
}

void Matrix::set_row(const int& row, const matrix_tuple& t) {
    gsl_matrix_set(__matrix->first, row, 0, std::get<0>(t));
    gsl_matrix_set(__matrix->first, row, 1, std::get<1>(t));
    gsl_matrix_set(__matrix->first, row, 2, std::get<2>(t));
    gsl_vector_set(__matrix->second, row, std::get<3>(t));
}

void Matrix::rotate(Coordinate& coord) const {
    gsl_vector* vec1 = gsl_vector_alloc(3);
    gsl_vector* vec2 = gsl_vector_alloc(3);
    gsl_vector_set(vec1, 0, coord.x());
    gsl_vector_set(vec1, 1, coord.y());
    gsl_vector_set(vec1, 2, coord.z());
    gsl_blas_dgemv(CblasNoTrans, 1, __matrix->first, vec1, 0, vec2);
    gsl_vector_add(vec2, __matrix->second);
    coord.set_x(gsl_vector_get(vec2, 0));
    coord.set_y(gsl_vector_get(vec2, 1));
    coord.set_z(gsl_vector_get(vec2, 2));
    gsl_vector_free(vec1);
    gsl_vector_free(vec2);
}

void Matrix::inverse_rotate(Coordinate& coord) const {
    gsl_vector* vec1 = gsl_vector_alloc(3);
    gsl_vector* vec2 = gsl_vector_alloc(3);
    gsl_vector_set(vec1, 0, coord.x());
    gsl_vector_set(vec1, 1, coord.y());
    gsl_vector_set(vec1, 2, coord.z());
    gsl_vector_sub(vec1, __matrix->second);
    gsl_blas_dgemv(CblasTrans, 1, __matrix->first, vec1, 0, vec2);
    coord.set_x(gsl_vector_get(vec2, 0));
    coord.set_y(gsl_vector_get(vec2, 1));
    coord.set_z(gsl_vector_get(vec2, 2));
    gsl_vector_free(vec1);
    gsl_vector_free(vec2);
}

std::ostream& operator<<(std::ostream& stream, const Matrix& m) {
    gsl_matrix* r = m.__matrix->first;
    gsl_vector* v = m.__matrix->second;
    stream << "[" << gsl_matrix_get(r, 0, 0) << ","
           << gsl_matrix_get(r, 0, 1) << "," << gsl_matrix_get(r, 0, 2)
           << "," << gsl_vector_get(v, 0) << "]" << std::endl;
    stream << "[" << gsl_matrix_get(r, 1, 0) << ","
           << gsl_matrix_get(r, 1, 1) << "," << gsl_matrix_get(r, 1, 2)
           << "," << gsl_vector_get(v, 1) << "]" << std::endl;
    stream << "[" << gsl_matrix_get(r, 2, 0) << ","
           << gsl_matrix_get(r, 2, 1) << "," << gsl_matrix_get(r, 2, 2)
           << "," << gsl_vector_get(v, 2) << "]" << std::endl;
    return stream;
}

}
}
