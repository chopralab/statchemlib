/* This is residue.hpp and is part of StatChemLIB
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

#ifndef RESIDUE_H
#define RESIDUE_H
#include "statchem/geometry/geometry.hpp"
#include "statchem/molib/atom.hpp"
#include "statchem/molib/it.hpp"

namespace statchem {

namespace molib {
class Atom;
class Chain;

class Residue : public template_map_container<Atom, Residue, Chain, int> {
   public:
    typedef enum {
        notassigned = 1,
        protein = 2,
        nucleic = 4,
        modres = 8,
        ion = 16,
        water = 32,
        hetero = 64
    } res_type;

    typedef std::tuple<char, std::string, int, char> res_tuple2;
    typedef std::pair<int, char> res_pair;
    typedef std::vector<Residue*> Vec;
    typedef std::set<Residue*> Set;

   private:
    std::string __resn;
    int __resi;
    char __ins_code;
    Residue::res_type __rest;
    geometry::Coordinate __crd;  // geometric center
   public:
    Residue(std::string resn, int resi, char ins_code, res_type(rest))
        : __resn(resn), __resi(resi), __ins_code(ins_code), __rest(rest) {}
    Residue(const Residue& rhs)
        : __resn(rhs.__resn),
          __resi(rhs.__resi),
          __ins_code(rhs.__ins_code),
          __rest(rhs.__rest),
          __crd(rhs.__crd) {
        for (auto& atom : rhs) {
            dbgmsg("Copy constructor : residue");
            add(new Atom(atom));
        }
    }
    void init_bio(Residue& residue_asym, const geometry::Matrix& bio_rota);
    void rotate(const geometry::Matrix& rota, const bool inverse = false);

    geometry::Coordinate& crd() { return __crd; }
    const geometry::Coordinate& crd() const { return __crd; }
    void set_crd() {
        for (auto& atom : *this) {
            __crd = __crd + atom.crd();
        }
        __crd = __crd / this->size();
    }  // calculate geom center
    Atom& add(Atom* a) { return this->aadd(a->atom_number(), a, this); }
    std::string resn() const { return __resn; }
    void set_resn(const std::string& resn) { __resn = resn; }
    void set_resi(int resi) { __resi = resi; }
    int resi() const { return __resi; }
    char ins_code() const { return __ins_code; }
    Residue::res_type rest() const { return __rest; }
    Atom& atom(int p) const { return this->element(p); }
    bool has_atom(int p) { return this->has_element(p); }
    Atom::Vec get_atoms(bool include_ions = true) const;
    void renumber_atoms(int new_start);
    Residue& erase_properties() {
        for (auto& atom : *this) atom.erase_properties();
        return *this;
    }

    Residue& regenerate_bonds(const Residue&);

    // NOTE: implementation in hydrogens.cpp
    void compute_hydrogen();
    void erase_hydrogen(bool temp_only = false);

    friend std::ostream& operator<<(std::ostream& stream, const Residue& r);
};

}  // molib
}

#endif
