/* This is capi.cpp and is part of StatChemLIB
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

#include "statchem/capi/capi.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/modeler/forcefield.hpp"
#include "statchem/modeler/modeler.hpp"
#include "statchem/molib/molecule.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
#include "statchem/fileio/fileout.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/score/kbff.hpp"

#include <boost/filesystem/path.hpp>
#include <memory>
#include <regex>

using namespace statchem;

std::unique_ptr<statchem::molib::Molecules> __receptor;
std::unique_ptr<statchem::molib::Molecules> __ligand;
std::unique_ptr<statchem::score::KBFF> __score;
std::unique_ptr<statchem::molib::Atom::Grid> __gridrec;
std::unique_ptr<statchem::OMMIface::ForceField> __ffield;
std::unique_ptr<statchem::OMMIface::Modeler> __modeler;
std::string __error_string = "";

const char* cd_get_error() { return __error_string.c_str(); }

size_t initialize_complex(const char* filename) {
    if (__receptor != nullptr) {
        __receptor.reset();
    }

    if (__ligand != nullptr) {
        __ligand.reset();
    }

    try {
        parser::FileParser rpdb(filename, parser::protein_poses_only, 1);

        __receptor = std::unique_ptr<statchem::molib::Molecules>(
            new statchem::molib::Molecules);

        rpdb.parse_molecule(*__receptor);

        __receptor->compute_idatm_type();

        __gridrec = std::unique_ptr<statchem::molib::Atom::Grid>(
            new statchem::molib::Atom::Grid(__receptor->get_atoms()));

        parser::FileParser lpdb(filename, parser::docked_poses_only, 1);

        __ligand = std::unique_ptr<statchem::molib::Molecules>(
            new statchem::molib::Molecules);

        lpdb.parse_molecule(*__ligand);

        __ligand->compute_idatm_type();

        return 1;
    } catch (std::exception& e) {
        __error_string = std::string("Error in loading complex: ") + e.what();
        return 0;
    }
}

size_t write_complex(const char* filename) {
    if (__ligand == nullptr || __receptor == nullptr) {
        __error_string = std::string(
            "You must run initialize_ligand and initialize_receptor before writing a complex");
        return 0;
    }

    try {

        double score = 0;
        if (__score != nullptr) {
            score =  calculate_score();
        }

        // output docked molecule conformations
        std::stringstream ss;
        fileio::print_complex_pdb(
            ss, (*__ligand)[0], (*__receptor)[0], score);
        fileio::output_file(ss.str(), filename);

        return 1;
    } catch (std::exception& e) {
        __error_string = std::string("Error in saving complex: ") + e.what();
        return 0;
    }
}

size_t initialize_receptor(const char* filename) {
    if (__receptor != nullptr) {
        __receptor.reset();
    }

    try {
        statchem::parser::FileParser rpdb(filename, statchem::parser::first_model,
                                         1);

        __receptor = std::unique_ptr<statchem::molib::Molecules>(
            new statchem::molib::Molecules);

        rpdb.parse_molecule(*__receptor);

        __receptor->compute_idatm_type()
            .compute_hydrogen()
            .compute_bond_order()
            .compute_bond_gaff_type()
            .refine_idatm_type()
            .erase_hydrogen()    // needed because refine changes connectivities
            .compute_hydrogen()  // needed because refine changes connectivities
            .compute_ring_type()
            .compute_gaff_type()
            .compute_rotatable_bonds()  // relies on hydrogens being assigned
            .erase_hydrogen();

        __gridrec = std::unique_ptr<statchem::molib::Atom::Grid>(
            new statchem::molib::Atom::Grid(__receptor->get_atoms()));

        return 1;
    } catch (std::exception& e) {
        __error_string = std::string("Error in loading receptor ") + e.what();
        return 0;
    }
}

size_t receptor_atom_count() {
    if (__receptor == nullptr) {
        __error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    return __receptor->element(0).get_atoms().size();
}

size_t receptor_atoms(size_t* idx, float* pos) {
    if (__receptor == nullptr) {
        __error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {
        statchem::molib::Atom::Vec atoms = __receptor->element(0).get_atoms();
        for (size_t i = 0; i < atoms.size(); ++i) {
            idx[i] = atoms[i]->atom_number();

            geometry::Point& crd = atoms[i]->crd();
            pos[i * 3 + 0] = crd.x();
            pos[i * 3 + 1] = crd.y();
            pos[i * 3 + 2] = crd.z();
        }

        return atoms.size();
    } catch (std::exception& e) {
        __error_string =
            std::string("Error creating receptor atom arrays: ") + e.what();
        return 0;
    }
}

size_t receptor_atom_details(char* chain_ids, size_t* resi, size_t* rest,
                             char* resn, size_t* elements) {
    if (__receptor == nullptr) {
        __error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {
        statchem::molib::Atom::Vec atoms = __receptor->element(0).get_atoms();
        for (size_t i = 0; i < atoms.size(); ++i) {
            chain_ids[i] = atoms[i]->br().br().chain_id();
            resi[i] = atoms[i]->br().resi();
            rest[i] = atoms[i]->br().rest();

            if (help::one_letter.find(atoms[i]->br().resn()) !=
                help::one_letter.end()) {
                resn[i] = help::one_letter.at(atoms[i]->br().resn());
            } else {
                resn[i] = 'X';
            }

            elements[i] = atoms[i]->element().number();
        }

        return atoms.size();
    } catch (std::exception& e) {
        __error_string =
            std::string("Error creating receptor atom arrays: ") + e.what();
        return 0;
    }
}

size_t receptor_bond_count() {
    if (__receptor == nullptr) {
        __error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {
        statchem::molib::BondSet bondset =
            statchem::molib::get_bonds_in(__receptor->get_atoms());
        return bondset.size();
    } catch (std::exception& e) {
        __error_string =
            std::string("Error creating receptor bond arrays: ") + e.what();
        return 0;
    }
}

size_t receptor_bonds(size_t* bonds) {
    if (__receptor == nullptr) {
        __error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {
        statchem::molib::BondSet bondset =
            statchem::molib::get_bonds_in(__receptor->get_atoms());

        size_t i = 0;

        for (const auto a : bondset) {
            bonds[i * 3 + 0] = a->atom1().atom_number();
            bonds[i * 3 + 1] = a->atom2().atom_number();

            bonds[i * 3 + 2] = 0;
            bonds[i * 3 + 2] |= a->is_single() ? SINGLE_BOND : 0;
            bonds[i * 3 + 2] |= a->is_double() ? DOUBLE_BOND : 0;
            bonds[i * 3 + 2] |= a->is_triple() ? TRIPLE_BOND : 0;
            bonds[i * 3 + 2] |= a->is_ring() ? INRING_BOND : 0;
            bonds[i * 3 + 2] |= a->is_rotatable() ? ROTATE_BOND : 0;
            bonds[i * 3 + 2] |= a->is_aromatic() ? AROMAT_BOND : 0;

            ++i;
        }

        return i;
    } catch (std::exception& e) {
        __error_string =
            std::string("Error creating receptor bond arrays: ") + e.what();
        return 0;
    }
}

size_t initialize_ligand(const char* filename) {
    if (__ligand != nullptr) {
        __ligand.reset();
    }

    try {
        statchem::parser::FileParser rpdb(filename, statchem::parser::first_model,
                                         1);

        __ligand = std::unique_ptr<statchem::molib::Molecules>(
            new statchem::molib::Molecules);

        rpdb.parse_molecule(*__ligand);

        __ligand->compute_idatm_type()
            .compute_hydrogen()
            .compute_bond_order()
            .compute_bond_gaff_type()
            .refine_idatm_type()
            .erase_hydrogen()    // needed because refine changes connectivities
            .compute_hydrogen()  // needed because refine changes connectivities
            .compute_ring_type()
            .compute_gaff_type()
            .compute_rotatable_bonds()  // relies on hydrogens being assigned
            .erase_hydrogen();

        return 1;
    } catch (std::exception& e) {
        __error_string = std::string("Error in loading ligand ") + e.what();
        return 0;
    }
}

size_t ligand_atom_count() {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    return __ligand->element(0).get_atoms().size();
}

size_t ligand_atoms(size_t* idx, float* pos) {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        statchem::molib::Atom::Vec atoms = __ligand->element(0).get_atoms();
        for (size_t i = 0; i < atoms.size(); ++i) {
            idx[i] = atoms[i]->atom_number();

            geometry::Point& crd = atoms[i]->crd();
            pos[i * 3 + 0] = crd.x();
            pos[i * 3 + 1] = crd.y();
            pos[i * 3 + 2] = crd.z();
        }

        return atoms.size();
    } catch (std::exception& e) {
        __error_string = std::string("Error creating atom arrays");
        return 0;
    }
}

size_t ligand_atom_details(char* chain_ids, size_t* resi, size_t* rest,
                           size_t* elements, int* idatm) {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        statchem::molib::Atom::Vec atoms = __ligand->element(0).get_atoms();
        for (size_t i = 0; i < atoms.size(); ++i) {
            chain_ids[i] = atoms[i]->br().br().chain_id();
            resi[i] = atoms[i]->br().resi();
            rest[i] = atoms[i]->br().rest();
            elements[i] = atoms[i]->element().number();
            idatm[i] = atoms[i]->idatm_type();
        }

        return atoms.size();
    } catch (std::exception& e) {
        __error_string =
            std::string("Error creating ligand atom arrays: ") + e.what();
        return 0;
    }
}

size_t ligand_ligand_atoms() {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
	}

	auto all_atoms = __ligand->get_atoms();
    __ligand->erase_properties();
    for (auto atom : all_atoms) {
        atom->set_idatm_type("???");
	}

	__ligand->compute_idatm_type();
    return 1;
}

size_t ligand_bond_count() {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        statchem::molib::BondSet bondset =
            statchem::molib::get_bonds_in(__ligand->get_atoms());
        return bondset.size();
    } catch (std::exception& e) {
        __error_string =
            std::string("Error creating ligand bond arrays: ") + e.what();
        return 0;
    }
}

size_t ligand_bonds(size_t* bonds) {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        statchem::molib::BondSet bondset =
            statchem::molib::get_bonds_in(__ligand->get_atoms());

        size_t i = 0;

        for (const auto a : bondset) {
            bonds[i * 3 + 0] = a->atom1().atom_number();
            bonds[i * 3 + 1] = a->atom2().atom_number();

            bonds[i * 3 + 2] = 0;
            bonds[i * 3 + 2] |= a->is_single() ? SINGLE_BOND : 0;
            bonds[i * 3 + 2] |= a->is_double() ? DOUBLE_BOND : 0;
            bonds[i * 3 + 2] |= a->is_triple() ? TRIPLE_BOND : 0;
            bonds[i * 3 + 2] |= a->is_ring() ? INRING_BOND : 0;
            bonds[i * 3 + 2] |= a->is_rotatable() ? ROTATE_BOND : 0;

            bonds[i * 3 + 2] |= a->is_aromatic() ? AROMAT_BOND : 0;

            ++i;
        }

        return i;
    } catch (std::exception& e) {
        __error_string =
            std::string("Error creating ligand bond arrays: ") + e.what();
        return 0;
    }
}

size_t ligand_get_neighbors(size_t atom_idx, size_t* neighbors) {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        statchem::molib::Atom::Vec atoms = __ligand->element(0).get_atoms();

        if (atom_idx >= atoms.size()) {
            __error_string = std::string("Atom index out of bounds");
            return 0;
        }

        const statchem::molib::Atom* current_atom = atoms[atom_idx];

        statchem::molib::BondVec bdev = current_atom->get_bonds();

        for (size_t i = 0; i < bdev.size(); ++i) {
            const auto& neigh = bdev[i]->second_atom(*current_atom);
            neighbors[i] = distance(atoms.begin(),
                                    find(atoms.begin(), atoms.end(), &neigh));
        }

        return bdev.size();
    } catch (std::exception& e) {
        __error_string = std::string("Error creating atom arrays");
        return 0;
    }
}

size_t initialize_scoring(const char* obj_dir) {
    return initialize_scoring_full(obj_dir, "radial", "mean", "complete", 15.0,
                                   0.01, 10.0);
}

size_t initialize_scoring_full(const char* obj_dir, const char* ref,
                               const char* func, const char* comp, float cutoff,
                               float step, float) {
    if (__ligand == nullptr || __receptor == nullptr) {
        __error_string = std::string(
            "You must run initialize_ligand and initialize_receptor first");
        return 0;
    }

    try {
        __score = std::unique_ptr<statchem::score::KBFF>(
            new statchem::score::KBFF(func, comp, ref, cutoff, step));

        boost::filesystem::path p(obj_dir);
        p /= "csd_complete_distance_distributions.txt.xz";

        __score
            ->define_composition(__receptor->get_idatm_types(),
                                 __ligand->get_idatm_types())
            .process_distributions(p.string())
            .compile_scoring_function();
        __score->compile_objective_function();

        return 1;
    } catch (std::exception& e) {
        __error_string =
            std::string("Error in creating scoring function: ") + e.what();
        return 0;
    }
}

size_t initialize_plugins(const char* plugin_dir) {
    try {
        statchem::OMMIface::SystemTopology::loadPlugins(plugin_dir);

        return 1;
    } catch (std::exception& e) {
        __error_string = std::string("Error in loading plugins: ") + e.what();
        return 0;
    }
}

size_t initialize_ffield(const char* data_dir, double dist_cutoff) {
    if (__score == nullptr) {
        __error_string = std::string("You must run initialize_score first");
        return 0;
    }

    try {
        boost::filesystem::path p(data_dir);

        __ffield = std::unique_ptr<statchem::OMMIface::ForceField>(
            new statchem::OMMIface::ForceField);

        __ffield->parse_gaff_dat_file((p / "gaff.dat").string())
            .add_kb_forcefield(*__score, dist_cutoff)
            .parse_forcefield_file((p / "amber10.xml").string())
            .parse_forcefield_file((p / "tip3p.xml").string());

        __receptor->element(0).prepare_for_mm(*__ffield, *__gridrec);

        __ffield->insert_topology(__receptor->element(0));
        __ffield->insert_topology(__ligand->element(0));

        return 1;
    } catch (std::exception& e) {
        __error_string =
            std::string("Error in creating forcefield: ") + e.what();
        return 0;
    }
}

size_t initialize_modeler(const char* platform, const char* precision,
                          const char* accelerators) {
    if (__ffield == nullptr) {
        __error_string = std::string("You must run initialize_ffield first");
        return 0;
    }
    try {
        __modeler = std::unique_ptr<OMMIface::Modeler>(
            new OMMIface::Modeler(*__ffield));

        auto rec_atoms = __receptor->get_atoms();
        auto lig_atoms = __ligand->get_atoms();

        __modeler->add_topology(rec_atoms);
        __modeler->add_topology(lig_atoms);

        __modeler->init_openmm(platform, precision, accelerators);

        return 1;
    } catch (std::exception& e) {
        __error_string = std::string("Error in creating modeler: ") + e.what();
        return 0;
    }
}

float calculate_score() {
    if (__score == nullptr) {
        __error_string = std::string("You must run initialize_score first");
        return 0;
    }

    try {
        return __score->non_bonded_energy(*__gridrec, (*__ligand)[0]);
    } catch (std::exception& e) {
        __error_string = std::string("Error in scoring: ") + e.what();
        return 0;
    }
}

size_t set_positions_ligand(const size_t* atoms, const float* positions,
                            size_t size) {
    if (__ligand == nullptr) {
        __error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        statchem::molib::Residue* residue =
            __ligand->element(0).get_residues().at(0);

        for (size_t i = 0; i < size; ++i) {
            residue->element(atoms[i]).set_crd(
                geometry::Point(positions[i * 3 + 0], positions[i * 3 + 1],
                                positions[i * 3 + 2]));
        }

        return 1;
    } catch (std::exception& e) {
        __error_string =
            std::string("Error in setting ligand coordinates: ") + e.what();
        return 0;
    }
}

size_t set_positions_receptor(const size_t* atoms, const float* positions,
                              size_t size) {
    if (__receptor == nullptr) {
        __error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {
        size_t resi = 0;
        statchem::molib::Residue* residue =
            __receptor->element(0).get_residues().at(resi++);

        for (size_t i = 0; i < size; ++i) {
            while (!residue->has_element(atoms[i])) {
                residue = __receptor->element(0).get_residues().at(resi++);
            }

            if (!residue->has_element(atoms[i])) {
                __error_string = std::string("Error: could not atom = ") +
                                 std::to_string(atoms[i]);
                break;
            }

            residue->element(atoms[i]).set_crd(
                geometry::Point(positions[i * 3 + 0], positions[i * 3 + 1],
                                positions[i * 3 + 2]));
        }

        __gridrec.reset(
            new statchem::molib::Atom::Grid(__receptor->get_atoms()));

        return 1;
    } catch (std::exception& e) {
        __error_string =
            std::string("Error in setting receptor coordinates: ") + e.what();
        return 0;
    }
}

size_t is_adjacent(size_t atom1, size_t atom2) {
    auto& residue = __ligand->first().first().first().first().first();

    auto& atom1_obj = residue.atom(atom1);
    auto& atom2_obj = residue.atom(atom2);

    return atom1_obj.is_adjacent(atom2_obj) || atom2_obj.is_adjacent(atom1_obj);
}

size_t add_ligand_bond(size_t atom1, size_t atom2) {
    auto& residue = __ligand->first().first().first().first().first();

    auto& atom1_obj = residue.atom(atom1);
    auto& atom2_obj = residue.atom(atom2);

    switch (atom1_obj.atom_number()) {
        case 6:
        case 7:
            if (atom1_obj.size() >= 4) {
                return 0;
            }
            break;
        case 8:
            if (atom1_obj.size() >= 2) {
                return 0;
            }
            break;
        case 9:
        case 17:
        case 35:
        case 53:
            if (atom1_obj.size() >= 1) {
                return 0;
            }
            break;
        default:
            break;
    }

    switch (atom2_obj.atom_number()) {
        case 6:
        case 7:
            if (atom2_obj.size() >= 4) {
                return 0;
            }
            break;
            if (atom2_obj.size() >= 4) {
                return 0;
            }
            break;
        case 8:
            if (atom2_obj.size() >= 2) {
                return 0;
            }
            break;
        case 9:
        case 17:
        case 35:
        case 53:
            if (atom2_obj.size() >= 1) {
                return 0;
            }
            break;
        default:
            break;
    }

    atom1_obj.connect(atom2_obj);

    return atom1_obj.is_adjacent(atom2_obj);
}

size_t remove_ligand_bond(size_t atom1, size_t atom2) {
    auto& residue = __ligand->first().first().first().first().first();

    auto& atom1_obj = residue.atom(atom1);
    auto& atom2_obj = residue.atom(atom2);

    atom1_obj.erase_bond(atom2_obj);
    atom2_obj.erase_bond(atom1_obj);

    size_t atom_to_remove = -1;
    for (size_t i = 0; i < atom1_obj.size(); ++i) {
        if (atom2_obj.atom_number() == atom1_obj[i].atom_number()) {
            atom_to_remove = i;
            break;
        }
    }

    if (atom_to_remove == (size_t)-1) {
        return 0;
    }

    atom1_obj.erase(atom_to_remove);

    for (size_t i = 0; i < atom2_obj.size(); ++i) {
        if (atom1_obj.atom_number() == atom2_obj[i].atom_number()) {
            atom_to_remove = i;
            break;
        }
    }

    atom2_obj.erase(atom_to_remove);

    return (!atom1_obj.is_adjacent(atom2_obj) && !atom2_obj.is_adjacent(atom1_obj));
}

size_t minimize_complex(size_t max_iter) {
    if (__modeler == nullptr) {
        __error_string = std::string("You must run initialize_modeler first");
        return 0;
    }

    try {
        __modeler->set_max_iterations(max_iter);

        auto rec_atoms = __receptor->get_atoms();
        auto lig_atoms = __ligand->get_atoms();

        __modeler->add_crds(rec_atoms, __receptor->get_crds());
        __modeler->add_crds(lig_atoms, __ligand->get_crds());

        __modeler->init_openmm_positions();

        __modeler->minimize_state();

        // init with minimized coordinates
        geometry::Point::Vec rec_coords = __modeler->get_state(rec_atoms);

        for (size_t i = 0; i < rec_atoms.size(); ++i)
            rec_atoms[i]->set_crd(rec_coords[i]);

        geometry::Point::Vec lig_coords = __modeler->get_state(lig_atoms);

        for (size_t j = 0; j < lig_atoms.size(); ++j)
            lig_atoms[j]->set_crd(lig_coords[j]);

        __gridrec.reset(
            new statchem::molib::Atom::Grid(__receptor->get_atoms()));

        return 1;
    } catch (std::exception& e) {
        __error_string =
            std::string("Error in setting receptor coordinates: ") + e.what();
        return 0;
    }
}
