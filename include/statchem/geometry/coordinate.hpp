/* This is coordinate.hpp and is part of StatChemLIB
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

#ifndef COORDINATE_H
#define COORDINATE_H
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "statchem/helper/error.hpp"

namespace statchem {

namespace geometry {
class Matrix;
class Coordinate {
    double __x, __y, __z;

   public:
    typedef std::vector<Coordinate> Vec;
    typedef std::set<Coordinate*> Set;
    typedef std::set<const Coordinate*> ConstSet;

    Coordinate() : __x(0), __y(0), __z(0) {}
    Coordinate(double x, double y, double z) : __x(x), __y(y), __z(z) {}
    Coordinate(const double* array)
        : __x(array[0]), __y(array[1]), __z(array[2]) {}
    Coordinate(const Coordinate& other)
        : __x(other.__x), __y(other.__y), __z(other.__z) {}

    void set_x(const double& x) { this->__x = x; }
    void set_y(const double& y) { this->__y = y; }
    void set_z(const double& z) { this->__z = z; }
    double x() const { return __x; }
    double y() const { return __y; }
    double z() const { return __z; }
    int i() const { return (int)::floor(__x); }
    int j() const { return (int)::floor(__y); }
    int k() const { return (int)::floor(__z); }
    const Coordinate& crd() const { return *this; }
    Coordinate floor() {
        return Coordinate(::floor(__x), ::floor(__y), ::floor(__z));
    }
    bool operator==(const Coordinate& right) const {
        return (fabs(this->x() - right.x()) <
                    std::numeric_limits<float>::epsilon() &&
                fabs(this->y() - right.y()) <
                    std::numeric_limits<float>::epsilon() &&
                fabs(this->z() - right.z()) <
                    std::numeric_limits<float>::epsilon());
    }
    bool operator<(const Coordinate& right) const {
        return (!(*this == right)) &&
               std::tie(this->__x, this->__y, this->__z) <
                   std::tie(right.__x, right.__y, right.__z);
    }
    bool operator>(const Coordinate& right) const {
        return (!(*this == right)) &&
               std::tie(this->__x, this->__y, this->__z) >
                   std::tie(right.__x, right.__y, right.__z);
    }
    static Coordinate cross(const Coordinate& l, const Coordinate& r) {
        return Coordinate(l.y() * r.z() - l.z() * r.y(),
                          -l.x() * r.z() + l.z() * r.x(),
                          l.x() * r.y() - l.y() * r.x());
    }
    static double scalar(const Coordinate& l, const Coordinate& r) {
        return l.x() * r.x() + l.y() * r.y() + l.z() * r.z();
    }
    Coordinate operator+(const double& right) const {
        return Coordinate(__x + right, __y + right, __z + right);
    }
    Coordinate operator+(const Coordinate& right) const {
        return Coordinate(__x + right.x(), __y + right.y(), __z + right.z());
    }
    Coordinate operator-(const double& right) const {
        return Coordinate(__x - right, __y - right, __z - right);
    }
    Coordinate operator-(const Coordinate& right) const {
        return Coordinate(__x - right.x(), __y - right.y(), __z - right.z());
    }
    Coordinate operator/(const double& right) const {
        if (right == 0)
            throw Error("Coordinate::operator/  division by zero\n");
        return Coordinate(__x / right, __y / right, __z / right);
    }
    Coordinate operator*(const double& right) const {
        return Coordinate(__x * right, __y * right, __z * right);
    }
    Coordinate operator-() const { return Coordinate(-x(), -y(), -z()); }
    double distance(const Coordinate& c) const {
        return sqrt(pow(__x - c.x(), 2) + pow(__y - c.y(), 2) +
                    pow(__z - c.z(), 2));
    }
    double distance_sq(const Coordinate& c) const {
        return pow(__x - c.x(), 2) + pow(__y - c.y(), 2) + pow(__z - c.z(), 2);
    }

    void distance(double) const {}  // just dummy : needed by grid

    void normalize() {
        const double length = this->distance(Coordinate(0.0, 0.0, 0.0));
        if (length == 0)
            throw Error("Coordinate::operator/  division by zero\n");
        this->__x /= length;
        this->__y /= length;
        this->__z /= length;
    }
    Coordinate norm() {
        Coordinate c = *this;
        c.normalize();
        return c;
    }
    std::string pdb() const {
        std::stringstream outs;
        outs << std::fixed << std::right << std::setprecision(3) << std::setw(8)
             << __x << std::fixed << std::right << std::setprecision(3)
             << std::setw(8) << __y << std::fixed << std::right
             << std::setprecision(3) << std::setw(8) << __z;
        return outs.str();
    }
    std::string simple() const {
        std::stringstream outs;
        outs << std::fixed << std::setprecision(3) << __x << " " << std::fixed
             << std::setprecision(3) << __y << " " << std::fixed
             << std::setprecision(3) << __z;
        return outs.str();
    }
    std::string with_underscores() const {
        std::stringstream outs;
        outs << std::setprecision(3) << __x << "_" << std::setprecision(3)
             << __y << "_" << std::setprecision(3) << __z;
        return outs.str();
    }
    std::unique_ptr<Coordinate> rotate(const Matrix&);
    std::unique_ptr<Coordinate> inverse_rotate(const Matrix&);
    friend std::ostream& operator<<(std::ostream& stream, const Coordinate& c) {
        stream << std::setprecision(8) << "[" << c.x() << "," << c.y() << ","
               << c.z() << "]";
        return stream;
    }
};
}
}

#endif
