/* This is atom.hpp and is part of StatChemLIB
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

#ifndef ATOM_H
#define ATOM_H
#include "statchem/geometry/coordinate.hpp"
#include "statchem/geometry/matrix.hpp"
#include "statchem/graph/graph.hpp"
#include "statchem/molib/bond.hpp"
#include "statchem/molib/element.hpp"
#include "statchem/molib/grid.hpp"
#include "statchem/molib/it.hpp"

namespace statchem {

namespace molib {
class Residue;

class Atom : public template_vector_container<Atom*, Atom> {
   public:
    typedef std::tuple<int, std::string, std::unique_ptr<geometry::Coordinate>,
                       double>
        atom_tuple;
    typedef std::pair<Atom*, Atom*> Pair;
    typedef std::vector<Atom*> Vec;
    typedef std::vector<const Atom*> ConstVec;
    typedef std::set<Atom*> Set;
    typedef std::set<const Atom*> CSet;
    typedef statchem::molib::Grid<Atom> Grid;
    typedef graph::Graph<Atom> Graph;

   private:
    int __atom_number;
    std::string __atom_name;
    geometry::Coordinate __crd;
    int __idatm_type;
    std::string __sybyl_type;
    std::string __gaff_type;
    Element __element;
    std::string __smiles_label;
    std::map<std::string, int> __smiles_prop;
    std::map<int, int> __aps;
    void* __br;  // back reference
    std::map<const Atom*, std::shared_ptr<Bond>> __bonds;

   public:
    Atom(const Atom& rhs)
        : __atom_number(rhs.__atom_number),
          __atom_name(rhs.__atom_name),
          __crd(rhs.__crd),
          __idatm_type(rhs.__idatm_type),
          __sybyl_type(rhs.__sybyl_type),
          __gaff_type(rhs.__gaff_type),
          __element(rhs.__element),
          __smiles_label(rhs.__smiles_label),
          __smiles_prop(rhs.__smiles_prop),
          __aps(rhs.__aps),
          __br(nullptr) {
        dbgmsg("Copy constructor : atom");
    }

    Atom(const int atom_number, const std::string& smiles_label,
         const std::map<std::string, int> smiles_prop)
        : __atom_number(atom_number),
          __element(""),
          __smiles_label(smiles_label),
          __smiles_prop(smiles_prop),
          __br(nullptr) {}
    Atom(const geometry::Coordinate& crd) : __crd(crd), __element("") {}
    Atom(const geometry::Coordinate& crd, const int& idatm_type)
        : __crd(crd), __idatm_type(idatm_type), __element("") {}
    // if element is missing, try to guess it from atom name
    Atom(int atom_number, const std::string& atom_name,
         const geometry::Coordinate& crd, const int idatm_type,
         const std::string& element = "", const std::string& sybyl_type = "")
        : __atom_number(atom_number),
          __atom_name(atom_name),
          __crd(crd),
          __idatm_type(idatm_type),
          __sybyl_type(sybyl_type),
          __gaff_type("???"),
          __element(element.empty() ? atom_name : element) {}
    Bond& connect(Atom& a2);
    void clear_bonds() { __bonds.clear(); }
    BondVec get_bonds() const {
        BondVec bonds;
        for (auto& kv : __bonds) bonds.push_back(&*kv.second);
        return bonds;
    }
    const std::shared_ptr<Bond>& get_shared_ptr_bond(const Atom& other) const {
        return __bonds.at(&other);
    }
    Bond& get_bond(const Atom& other) const {
        return *get_shared_ptr_bond(other);
    }
    std::shared_ptr<Bond>& insert_bond(const Atom& other,
                                       const std::shared_ptr<Bond>& bond) {
        return __bonds.insert({&other, bond}).first->second;
    }
    std::shared_ptr<Bond>& insert_bond(const Atom& other, Bond* bond) {
        return __bonds.insert({&other, std::shared_ptr<Bond>(bond)})
            .first->second;
    }
    void erase_bond(const Atom& other) { __bonds.erase(&other); }
    bool is_adjacent(const Atom& other) const {
        return __bonds.count(&other) != 0;
    }
    bool is_adjacent(const std::string& atom_name) const {
        for (auto& other : *this)
            if (other.atom_name() == atom_name) return true;
        return false;
    }
    int get_num_hydrogens() const;
    int atom_number() const { return __atom_number; }
    void set_atom_name(const std::string& atom_name) {
        __atom_name = atom_name;
    }
    void set_atom_number(int atom_number) { __atom_number = atom_number; }
    void set_idatm_type(const std::string& idatm_type);
    void set_gaff_type(const std::string& gaff_type) {
        __gaff_type = gaff_type;
    }
    Atom& insert_property(const std::string& prop, const int count) {
        __smiles_prop.insert({prop, count});
        return *this;
    }
    Atom& add_property(const std::string& prop) {
        __smiles_prop[prop]++;
        return *this;
    }
    Atom& erase_property(const std::string& prop) {
        __smiles_prop.erase(prop);
        return *this;
    }
    Atom& erase_properties() {
        __smiles_prop.clear();
        return *this;
    }
    bool has_property(const std::string& prop) const {
        return __smiles_prop.count(prop) != 0;
    }
    int get_num_property(const std::string& prop) const {
        dbgmsg("__smiles prop(" << prop
                                << ") count = " << __smiles_prop.count(prop));
        if (__smiles_prop.count(prop) == 0) return 0;
        return __smiles_prop.at(prop);
    }
    int get_num_bond_with_bond_gaff_type(const std::string& prop) const;
    int compute_num_property(const std::string& prop) const;
    void set_crd(const geometry::Coordinate& crd) { __crd = crd; }
    void rotate(const geometry::Matrix& matrix) { matrix.rotate(__crd); }
    void inverse_rotate(const geometry::Matrix& matrix) { matrix.inverse_rotate(__crd); }
    std::string idatm_type_unmask() const;
    int idatm_type() const { return __idatm_type; }
    const std::string& sybyl_type() const { return __sybyl_type; }
    double radius() const;
    const std::string& gaff_type() const { return __gaff_type; }
    const std::string& atom_name() const { return __atom_name; }
    Element element() const { return __element; }
    void set_element(const Element& e) { __element = e; };
    geometry::Coordinate& crd() { return __crd; }
    const geometry::Coordinate& crd() const { return __crd; }
    friend std::ostream& operator<<(std::ostream& stream, const Atom& a);
    double distance() const { return 0.0; }  // just dummy : needed by grid
    void distance(double) const {}           // just dummy : needed by grid
    const std::map<int, int>& get_aps() const { return __aps; }
    void set_members(const std::string& str);
    const Residue& br() const { return *static_cast<const Residue*>(__br); }
    void set_br(void* br) { __br = br; }
    // the following are required for Atom::Graph :-)
    bool compatible(const Atom& atom) const;
    std::string get_label() const;
    int weight() const { return 0; }  // dummy for graph ostream operator

    static Graph create_graph(const Vec& atoms);
    static Graph create_graph(const Set& atoms);
    static double compute_rmsd(const Graph& g1, const Graph& g2);
};

double compute_rmsd_vina_sq(const Atom::Vec& crds1, const Atom::Vec& crds2);
double compute_rmsd_vina(const Atom::Vec& crds1, const Atom::Vec& crds2);

BondSet get_bonds_in(const Atom::Set& atoms, bool in = true);
BondSet get_bonds_in(const Atom::Vec& atoms, bool in = true);

std::ostream& operator<<(std::ostream& stream, const Atom::Set& atoms);
std::ostream& operator<<(std::ostream& stream, const Atom::Vec& atoms);
}
}

#endif
