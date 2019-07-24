/* This is residue.cpp and is part of StatChemLIB
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

#include "statchem/molib/residue.hpp"
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "statchem/molib/atom.hpp"
#include "statchem/molib/atomtype.hpp"
#include "statchem/molib/bondtype.hpp"
using namespace std;

namespace statchem {
namespace molib {

ostream& operator<<(ostream& stream, const Residue& r) {
    for (auto& atom : r) {
        stream << atom;
    }
    return stream;
}

void Residue::init_bio(Residue& residue_asym,
                       const geometry::Matrix& bio_rota) {
    for (Atom& atom : residue_asym) {
        Atom& last = add(new Atom(atom));
        last.rotate(bio_rota);  // do the actual rotation
    }
}
void Residue::rotate(const geometry::Matrix& rota, const bool inverse) {
    if (inverse)
        for (auto& atom : *this) {
            atom.inverse_rotate(rota);
        }
    else
        for (auto& atom : *this) {
            atom.rotate(rota);
        }
}

Atom::Vec Residue::get_atoms(bool include_ions) const {
    Atom::Vec atoms;
    for (auto& atom : *this) {
        if (!include_ions && atom.element() >= Element::Sc &&
            atom.element() <= Element::Zn) {
            continue;
        }

        atoms.push_back(&atom);
    }
    return atoms;
}

void Residue::renumber_atoms(int new_start) {
    for (auto& atom : *this) {
        atom.set_atom_number(new_start++);
    }
}

Residue& Residue::regenerate_bonds(const Residue& template_model) {
    // copied molecule must be a fragment of template molecule...
    map<int, Atom*> atom_number_to_atom;
    for (auto& patom : template_model.get_atoms()) {
        atom_number_to_atom[patom->atom_number()] = patom;
    }
    map<const Atom*, Atom*> atom1_to_copy1;
    for (auto& patom : this->get_atoms()) {
        patom->clear();  // clear bondees as they refer to the template molecule
        patom->clear_bonds();  // clear bonds : fixes double CONECT words bug in
                               // PDB output
        atom1_to_copy1[atom_number_to_atom.at(patom->atom_number())] = patom;
    }
    for (auto& kv : atom1_to_copy1) {  // regenerate bonds
        const Atom& atom1 = *kv.first;
        Atom& copy1 = *kv.second;
        for (auto& atom2 : atom1) {
            if (atom1_to_copy1.count(&atom2)) {
                // copying of bonds does not preserve bond properties !!!
                Atom& copy2 = *atom1_to_copy1.at(&atom2);
                copy1.connect(copy2);
            }
        }
    }
    connect_bonds(get_bonds_in(this->get_atoms()));
    return *this;
}
}
}
