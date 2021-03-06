/* This is molecules.cpp and is part of StatChemLIB
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

#include "statchem/molib/molecule.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "statchem/fragmenter/fragmenter.hpp"
#include "statchem/fragmenter/unique.hpp"
#include "statchem/geometry/geometry.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/molib/assembly.hpp"
#include "statchem/molib/atom.hpp"
#include "statchem/molib/atomtype.hpp"
#include "statchem/molib/bond.hpp"
#include "statchem/molib/bondtype.hpp"
#include "statchem/molib/chain.hpp"
#include "statchem/molib/model.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/molib/residue.hpp"

using namespace std;

namespace statchem {
namespace molib {

ostream& operator<<(ostream& stream, const Molecules& m) {
    stream << setw(6) << left << "REMARK   4 NRPDB " << m.name() << endl;
    for (auto& molecule : m) {
        stream << molecule;
    }
    stream << setw(6) << left << "REMARK   4 END" << endl;
    return stream;
}

set<int> Molecules::get_idatm_types(set<int> previous) const {
    for (auto& molecule : *this)
        for (auto& pa : molecule.get_atoms()) {
            previous.insert(pa->idatm_type());
        }
    return previous;
}

void Molecules::rotate(const geometry::Matrix& rota, const bool inverse) {
    for (auto& molecule : *this) {
        molecule.rotate(rota, inverse);
    }
}

double Molecules::compute_max_radius() const {
    geometry::Coordinate gc = compute_geometric_center();
    double max_radius = -HUGE_VAL;
    for (auto& molecule : *this) {
        Atom::Vec atoms = molecule.get_atoms();
        for (auto& patom : atoms)
            max_radius = max(max_radius, patom->crd().distance(gc));
    }
    return max_radius;
}

geometry::Coordinate Molecules::compute_geometric_center() const {
    geometry::Coordinate center;
    for (auto& molecule : *this) {
        center = center + molecule.compute_geometric_center();
    }
    return center / this->size();
}

Molecules& Molecules::compute_idatm_type() {
    auto& mols = *this;
    for (size_t i = 0; i < mols.size(); ++i) {
        auto& molecule = mols[i];
        try {
            AtomType::compute_idatm_type(molecule.get_atoms());
        } catch (exception& e) {
            log_error << "errmesg : deleting molecule " << molecule.name()
                      << " (computing idatm types failed) due to " << e.what()
                      << endl;
            mols.erase_shrink(i--);
        }
    }
    return *this;
}

Molecules& Molecules::compute_hydrogen() {
    auto& mols = *this;
    for (size_t i = 0; i < mols.size(); ++i) {
        auto& molecule = mols[i];
        try {
            for (auto& presidue : molecule.get_residues()) {
                Residue& residue = *presidue;
                if (!(help::standard_residues.count(residue.resn()) ||
                      help::ions.count(residue.resn()))) {
                    residue.compute_hydrogen();
                }
            }
        } catch (exception& e) {
            log_error << "errmesg : deleting molecule " << molecule.name()
                      << " (computing hydrogens failed)" << endl;
            mols.erase_shrink(i--);
        }
    }
    return *this;
}

Molecules& Molecules::compute_bond_order() {
    auto& mols = *this;
    for (size_t i = 0; i < mols.size(); ++i) {
        auto& molecule = mols[i];
        try {
            for (auto& presidue : molecule.get_residues()) {
                Residue& residue = *presidue;
                if (!(help::standard_residues.count(residue.resn()) ||
                      help::ions.count(residue.resn()))) {
                    BondOrder::compute_bond_order(residue.get_atoms(false));
                }
            }
        } catch (exception& e) {
            log_error << "errmesg : deleting molecule " << molecule.name()
                      << " (bond order assignment failed) because " << e.what()
                      << endl;
            mols.erase_shrink(i--);
        }
    }
    return *this;
}

Molecules& Molecules::compute_bond_gaff_type() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        if (!(help::standard_residues.count(residue.resn()) ||
              help::ions.count(residue.resn()))) {
            BondOrder::compute_bond_gaff_type(residue.get_atoms(false));
        }
    }
    return *this;
}

Molecules& Molecules::refine_idatm_type() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        if (!(help::standard_residues.count(residue.resn()) ||
              help::cofactor_residues.count(residue.resn()) ||
              help::ions.count(residue.resn()))) {
            AtomType::refine_idatm_type(residue.get_atoms());
        }
    }
    return *this;
}

Molecules& Molecules::compute_chirality() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        if (!(help::standard_residues.count(residue.resn()) ||
              help::cofactor_residues.count(residue.resn()) ||
              help::ions.count(residue.resn()))) {
            statchem::molib::compute_chirality(residue.get_atoms());
        }
    }
    return *this;
}

Molecules& Molecules::erase_hydrogen() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        // Remove hydrogens added to protonated ions as well
        if (!(help::standard_residues.count(residue.resn()) ||
              help::ions.count(residue.resn()))) {
            residue.erase_hydrogen(false);
        }
    }
    return *this;
}

Molecules& Molecules::erase_temporary_hydrogen() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        // Remove hydrogens added to protonated ions as well
        if (!(help::standard_residues.count(residue.resn()) ||
              help::ions.count(residue.resn()))) {
            residue.erase_hydrogen(true);
        }
    }
    return *this;
}

Molecules& Molecules::compute_ring_type() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        if (!(help::standard_residues.count(residue.resn()) ||
              help::ions.count(residue.resn()))) {
            AtomType::compute_ring_type(residue.get_atoms());
        }
    }
    return *this;
}

Molecules& Molecules::compute_gaff_type() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        if (!(help::standard_residues.count(residue.resn()) ||
              help::metals.count(residue.resn()))) {
            AtomType::compute_gaff_type(residue.get_atoms(false));
            for (auto& atom : residue.get_atoms()) {
                if (atom->element() >= Element::Sc &&
                    atom->element() <= Element::Zn) {
                    std::string new_gaff_name = atom->element().name();
                    std::transform(new_gaff_name.begin(), new_gaff_name.end(),
                                   new_gaff_name.begin(), ::tolower);
                    atom->set_gaff_type(new_gaff_name);
                }

                if (residue.resn() == "HEM") {
                    if (atom->atom_name() == "NA" ||
                        atom->atom_name() == "NC") {
                        atom->set_gaff_type("nd");
                    }

                    if (atom->atom_name() == "NB" ||
                        atom->atom_name() == "ND") {
                        atom->set_gaff_type("nc");
                    }
                }
            }
        }
    }
    return *this;
}
Molecules& Molecules::compute_rotatable_bonds() {
    for (auto& presidue : this->get_residues()) {
        Residue& residue = *presidue;
        if (!(help::standard_residues.count(residue.resn()) ||
              help::cofactor_residues.count(residue.resn()) ||
              help::ions.count(residue.resn()))) {
            BondOrder::compute_rotatable_bonds(residue.get_atoms(false));
        }
    }
    return *this;
}
Molecule::Vec Molecules::get_molecules(const Residue::res_type& rest) const {
    Molecule::Vec molecules;
    for (auto& molecule : *this) {
        Residue& first = molecule.first().first().first().first();
        if (first.rest() == rest) molecules.push_back(&molecule);
    }
    return molecules;
}

Atom::Vec Molecules::get_atoms(const string& chain_ids,
                               const Residue::res_type& rest,
                               const int model_number) const {
    Atom::Vec atoms;
    for (auto& molecule : *this) {
        auto ret = molecule.get_atoms(chain_ids, rest, model_number);
        atoms.insert(atoms.end(), ret.begin(), ret.end());
    }
    return atoms;
}

Residue::Vec Molecules::get_residues() const {
    Residue::Vec residues;
    for (auto& molecule : *this)
        for (auto& assembly : molecule)
            for (auto& model : assembly)
                for (auto& chain : model)
                    for (auto& residue : chain) {
                        residues.push_back(&residue);
                    }
    return residues;
}

geometry::Point::Vec Molecules::get_crds(const string& chain_ids,
                                         const Residue::res_type& rest,
                                         const int model_number) const {
    geometry::Point::Vec crds;
    for (auto& patom : this->get_atoms(chain_ids, rest, model_number))
        crds.push_back(patom->crd());
    return crds;
}

Molecules& Molecules::compute_overlapping_rigid_segments(
    const string& seeds_file) {
    Unique u(seeds_file);
    for (auto& molecule : *this) molecule.compute_overlapping_rigid_segments(u);
    return *this;
}

void create_mols_from_seeds(set<int>& added, molib::Molecules& seeds,
                            const molib::Molecules& mols) {
    for (auto& molecule : mols)
        for (auto& assembly : molecule)
            for (auto& model : assembly) {
                for (auto& fragment :
                     model.get_rigid()) {  // iterate over seeds
                    if (fragment.is_seed()) {
                        dbgmsg("considering to add " << fragment.get_seed_id());
                        if (!added.count(
                                fragment.get_seed_id())) {  // take seeds that
                                                            // haven't been
                                                            // docked already
                            dbgmsg("added " << fragment.get_seed_id());
                            added.insert(fragment.get_seed_id());
                            // add to new molecules
                            molib::Molecule& seed =
                                seeds.add(new molib::Molecule(
                                    std::to_string(fragment.get_seed_id())));
                            molib::Assembly& a =
                                seed.add(new molib::Assembly(0));
                            molib::Model& mod = a.add(new molib::Model(1));
                            molib::Chain& c = mod.add(new molib::Chain('X'));
                            molib::Residue& r = c.add(new molib::Residue(
                                "XXX", 1, ' ', molib::Residue::hetero));
                            for (const molib::Atom* atom : fragment.get_all()) {
                                r.add(new molib::Atom(*atom));
                            }
                            seed.regenerate_bonds(molecule);
                        }
                    }
                }
            }
}
}
}
