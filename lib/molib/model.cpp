/* This is model.cpp and is part of StatChemLIB
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

#include "statchem/molib/model.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "statchem/molib/atom.hpp"
#include "statchem/molib/bond.hpp"
#include "statchem/molib/chain.hpp"
#include "statchem/molib/residue.hpp"
using namespace std;

namespace statchem {
namespace molib {
ostream& operator<<(ostream& stream, const Model& m) {
    stream << setw(6) << left << "MODEL" << setw(4) << " " << setw(4) << right
           << m.number() << endl;
    for (auto& r1 : m.__remarks) {
        const int& remark_number = r1.first.first;
        const Residue::res_tuple2& ligand = r1.first.second;
        stream << setw(6) << left << "REMARK" << setw(4) << right
               << remark_number << setw(6) << right << "ALRES"
               << " " << std::get<0>(ligand) << ":" << std::get<1>(ligand)
               << ":" << std::get<2>(ligand) << ":" << std::get<3>(ligand);
        for (auto& rpair : r1.second) {
            stream << "," << std::get<0>(rpair.first) << ":"
                   << std::get<1>(rpair.first) << ":"
                   << std::get<2>(rpair.first) << ":"
                   << std::get<3>(rpair.first);
            stream << "=" << std::get<0>(rpair.second) << ":"
                   << std::get<1>(rpair.second) << ":"
                   << std::get<2>(rpair.second) << ":"
                   << std::get<3>(rpair.second);
        }
        stream << endl;
    }
    for (auto& chain : m) {
        stream << chain;
    }
    for (auto& rigid : m.get_rigid()) {
        stream << "REMARK   8 RIGID " << rigid.get_seed_id();
        for (auto& atom : rigid.get_core())
            stream << " " << 'c' << atom->atom_number();
        for (auto& atom : rigid.get_join())
            stream << " " << 'j' << atom->atom_number();
        stream << endl;
    }
    for (auto& pbond : get_bonds_in(m.get_atoms())) {
        Bond& bond = *pbond;
        if (bond.is_rotatable()) {
            stream << "REMARK   8 ROTA " << bond.get_rotatable() << " "
                   << bond.atom1().atom_number() << " "
                   << bond.atom2().atom_number()
                   //<< " " << bond.stereo()
                   << endl;
        }
    }
    for (auto& pbond : get_bonds_in(m.get_atoms())) {
        Bond& bond = *pbond;
        stream << "REMARK   8 BONDTYPE " << bond.get_bo() << " "
               << bond.get_bond_gaff_type() << " " << bond.atom1().atom_number()
               << " " << bond.atom2().atom_number()
               //<< " " << bond.stereo()
               << endl;
    }
    stream << "ENDMDL" << endl;
    return stream;
}

void Model::init_bio(const Model& model_asym, const geometry::Matrix& matrix,
                     const set<char>& chains) {
    for (const char& chain_id : chains) {
        if (model_asym.has_chain(chain_id)) {
            Chain& last = add(new Chain(chain_id));
            last.init_bio(model_asym.chain(chain_id), matrix);
        }
    }
}
void Model::rotate(const geometry::Matrix& rota, const bool inverse) {
    for (auto& chain : *this) {
        chain.rotate(rota, inverse);
    }
}

Atom::Vec Model::get_atoms(const string& chain_ids,
                           const Residue::res_type& rest) const {
    Atom::Vec atoms;
    for (auto& chain : *this) {
        if (chain_ids == "" ||
            chain_ids.find(chain.chain_id()) != string::npos) {
            auto ret = chain.get_atoms(rest);
            atoms.insert(atoms.end(), ret.begin(), ret.end());
        }
    }
    return atoms;
}

Model& Model::regenerate_bonds(const Model& template_model) {
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
                Bond& new_bond = copy1.connect(copy2);
                Bond& old_bond = atom2.get_bond(atom1);
                new_bond.set_bo(old_bond.get_bo());
                new_bond.set_bond_gaff_type(old_bond.get_bond_gaff_type());
                // new_bond.set_stereo(old_bond.stereo());
            }
        }
    }
    connect_bonds(get_bonds_in(this->get_atoms()));
    return *this;
}
}
}
