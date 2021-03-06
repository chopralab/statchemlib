/* This is hydrogens.cpp and is part of StatChemLIB
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

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <regex>
#include <set>
#include <string>
#include <vector>
#include "statchem/fragmenter/fragmenter.hpp"
#include "statchem/fragmenter/unique.hpp"
#include "statchem/geometry/geometry.hpp"
#include "statchem/helper/benchmark.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/molib/molecule.hpp"
using namespace std;

namespace statchem {
namespace molib {
void Residue::compute_hydrogen() {
    try {
        if (this->empty()) {
            throw Error("die : compute hydrogens for empty atoms set?");
        }

        Atom::Vec all_atoms = this->get_atoms();

        // atom numbers must be unique for the WHOLE molecule, that is why
        // max atom number has to be found on molecule level
        int max_atom_number = 0;
        const Molecule& molecule = this->br().br().br().br();
        for (auto& patom : molecule.get_atoms())
            if (patom->atom_number() > max_atom_number)
                max_atom_number = patom->atom_number();

        for (auto& atom : *this) {
            dbgmsg("computing hydrogens for residue = " << this->resn());

            if (atom.element() !=
                Element::H) {  // don't visit "just added" hydrogens
                dbgmsg("hydro for : " << atom.idatm_type_unmask());
                int con =
                    help::get_info_map(atom.idatm_type_unmask()).substituents;

                /* EXCEPTION for 3-substituted P */
                if ((atom.idatm_type_unmask() == "P" ||
                     atom.idatm_type_unmask() == "PA") &&
                    atom.size() == 3)
                    con = 3;
                /* ***************************** */

                int num_h = con - atom.size();

                for (size_t i = 0; i < atom.size(); ++i)
                    if (atom[i].element() >= Element::Sc &&
                        atom[i].element() <= Element::Zn)
                        ++num_h;

                if (atom.element() >= Element::Sc &&
                    atom.element() <= Element::Zn)
                    num_h = 0;

                dbgmsg("computing hydrogens for "
                       << atom.idatm_type_unmask() << " con = " << con
                       << " atom.size() = " << atom.size());
                if (num_h > 0) {
                    dbgmsg("computing missing hydrogens for atom "
                           << atom << " con = " << con
                           << " atom.size() = " << atom.size()
                           << " number of added hydrogens = " << num_h);
                    // add dummy hydrogen atoms with zero coordinates
                    for (int i = 0; i < num_h; ++i) {
                        const string idatm_type =
                            (atom.element() == Element::C ? "HC" : "H");
                        dbgmsg("idatm_type = " << idatm_type);
                        dbgmsg(
                            "idatm_mask = " << help::idatm_mask.at(idatm_type));
                        Atom& hatom = this->add(new Atom(
                            ++max_atom_number, "H", geometry::Coordinate(),
                            help::idatm_mask.at(idatm_type)));
                        atom.connect(hatom);
                        dbgmsg("added hydrogen");
                        all_atoms.push_back(&hatom);
                    }
                } else if (num_h < 0) {
                    int h_excess = abs(num_h);
                    dbgmsg(
                        "deleting excess hydrogens because according to IDATM "
                        "type : "
                        << atom.idatm_type_unmask() << " this atom : " << atom
                        << " should have " << h_excess << " less hydrogens!");
                    // deleting hydrogens
                    for (size_t i = 0; i < atom.size(); ++i) {
                        auto& bondee = atom[i];
                        if (bondee.element() == Element::H) {
                            if (h_excess-- > 0) {
                                atom.erase(i--);
                                auto& shpbond =
                                    atom.get_shared_ptr_bond(bondee);
                                const Bond& deleted_bond = *shpbond;
                                erase_stale_bond_refs(
                                    deleted_bond,
                                    atom.get_bonds());  // delete references
                                dbgmsg(
                                    "shared_count1 = " << shpbond.use_count());
                                atom.erase_bond(bondee);
                                dbgmsg(
                                    "shared_count2 = " << shpbond.use_count());
                                bondee.erase_bond(atom);
                                dbgmsg(
                                    "shared_count3 = " << shpbond.use_count());
                                dbgmsg("residue before erasing hydrogen "
                                       << bondee.atom_number() << endl
                                       << *this);
                                this->erase(bondee.atom_number());
                                dbgmsg("residue after erasing hydrogen "
                                       << bondee.atom_number() << endl
                                       << *this);

                                auto it = find(all_atoms.begin(),
                                               all_atoms.end(), &bondee);
                                if (it != all_atoms.end()) {
                                    all_atoms.erase(it);
                                } else {
                                    throw Error(
                                        "die: cannot find newly added hydrogen "
                                        "in vector");
                                }
                            }
                        }
                    }
                    if (h_excess > 0) {
                        log_error << "Issue with " << '\n'
                                  << atom << "Excess of " << h_excess << endl;
                        throw Error(
                            "die : deleting of excess hydrogens failed");
                    }
                }
            }
        }
        // connect the new hydrogen bonds with the rest of the bond-graph
        erase_bonds(get_bonds_in(all_atoms));
        connect_bonds(get_bonds_in(all_atoms));
        dbgmsg("MOLECULE AFTER COMPUTING HYDROGENS " << endl << all_atoms);
        dbgmsg("BONDS AFTER COMPUTING HYDROGENS " << endl
                                                  << get_bonds_in(all_atoms));
    } catch (exception& e) {
        log_error << "errmesg : " << e.what() << endl;
        throw e;
    }
}
void Residue::erase_hydrogen(bool temp_only) {
    Atom::Vec natoms = this->get_atoms();
    natoms.erase(
        std::remove_if(
            natoms.begin(), natoms.end(),
            [this, temp_only](Atom* patom) -> bool {

                Atom& hatom = *patom;
                dbgmsg("deleting hydrogens for residue = " << this->resn());
                if (hatom.element() == Element::H &&
                     (!temp_only || hatom.crd() == statchem::geometry::Coordinate(0,0,0)) ) {
                    dbgmsg("deleting this hydrogen atom : " << hatom);
                    // H is always bonded to just one other heavy atom
                    auto& bondee = hatom[0];

                    // delete the H from bondee's neighbor list
                    for (size_t i = 0; i < bondee.size(); ++i) {
                        if (bondee[i].atom_number() == hatom.atom_number()) {
                            dbgmsg("erasing H from bondee at index = " << i);
                            bondee.erase(i);
                            break;
                        }
                    }
                    // erase bond to bondee-H
                    bondee.erase_bond(hatom);
                    // erase H from this residue
                    this->erase(hatom.atom_number());
                    return true;  // erase H from atoms vector
                }
                return false;

            }),
        natoms.end());
    dbgmsg("MOLECULE AFTER ERASING HYDROGENS " << endl << natoms);
}
}
}
