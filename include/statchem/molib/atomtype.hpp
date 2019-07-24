/* This is atomtype.hpp and is part of StatChemLIB
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

#ifndef ATOMTYPE_H
#define ATOMTYPE_H
#include <string>
#include <vector>
#include "statchem/molib/atom.hpp"
#include "statchem/molib/it.hpp"

namespace statchem {

namespace molib {
class Atom;
class Molecule;

namespace AtomType {
void compute_idatm_type(const Atom::Vec& atoms);
void refine_idatm_type(const Atom::Vec& atoms);
void compute_gaff_type(const Atom::Vec& atoms);
void compute_ring_type(const Atom::Vec& atoms);

std::tuple<double, size_t, size_t, size_t> determine_lipinski(
    const Atom::Vec& atoms);
}
}
}

#endif
