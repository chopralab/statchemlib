/* This is geometry.hpp and is part of StatChemLIB
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

#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <map>
#include <vector>
#include "statchem/geometry/coordinate.hpp"

namespace statchem {

namespace geometry {
typedef Coordinate Point;
typedef Coordinate Vector3;
typedef std::map<int, geometry::Point::Vec> GridPoints;

double degrees(double);
double radians(double);

double angle(const Vector3&, const Vector3&);
double angle(const Point&, const Point&, const Point&);
double dihedral(const Point&, const Point&, const Point&, const Point&);
Point line_evaluate(const Point&, const Vector3&, const double);

double compute_rmsd_sq(const Point::Vec& crds1, const Point::Vec& crds2);
double compute_rmsd(const Point::Vec& crds1, const Point::Vec& crds2);

Point compute_geometric_center(const geometry::Point::Vec& crds);
Point::Vec uniform_sphere(const int n);

std::ostream& operator<<(std::ostream& os, const geometry::Point::Vec& points);
std::ostream& operator<<(std::ostream& os,
                         const geometry::Point::ConstSet& points);
}
}

#endif
