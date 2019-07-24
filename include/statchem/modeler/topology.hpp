/* This is topology.hpp and is part of StatChemLIB
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

#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <map>
#include <set>
#include <string>
#include <vector>
#include "statchem/geometry/coordinate.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/molib/molecule.hpp"

namespace statchem {

namespace molib {
class Atom;
class Molecule;
}

namespace OMMIface {
struct ForceField;

class Topology {
   public:
    typedef std::set<std::pair<molib::Atom*, molib::Atom*>> BondedExclusions;
    typedef std::vector<std::pair<molib::Atom*, molib::Atom*>> Bonds;
    typedef std::vector<std::tuple<molib::Atom*, molib::Atom*, molib::Atom*>>
        Angles;
    typedef std::vector<
        std::tuple<molib::Atom*, molib::Atom*, molib::Atom*, molib::Atom*>>
        Dihedrals;
    molib::Atom::Vec atoms;
    Bonds bonds;
    Angles angles;
    Dihedrals dihedrals, impropers;
    BondedExclusions bonded_exclusions;

   private:
    std::map<const molib::Atom*, const int> atom_to_type;
    std::map<const molib::Atom*, const int> atom_to_index;

   public:
    ~Topology() { dbgmsg("calling destructor of Topology"); }

    Topology& add_topology(const molib::Atom::Vec& atoms,
                           const ForceField& ffield);
    int get_index(const molib::Atom& atom) const;
    int get_type(const molib::Atom& atom) const;
};

std::ostream& operator<<(std::ostream& stream,
                         const Topology::BondedExclusions& bonds);
std::ostream& operator<<(std::ostream& stream, const Topology::Bonds& bonds);
std::ostream& operator<<(std::ostream& stream, const Topology::Angles& angles);
std::ostream& operator<<(std::ostream& stream,
                         const Topology::Dihedrals& dihedrals);
}
}

#endif
