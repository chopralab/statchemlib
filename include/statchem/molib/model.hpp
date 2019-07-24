/* This is modeler.hpp and is part of StatChemLIB
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

#ifndef MODEL_H
#define MODEL_H
#include "statchem/fragmenter/fragmenter.hpp"
#include "statchem/geometry/geometry.hpp"
#include "statchem/geometry/matrix.hpp"
#include "statchem/molib/atom.hpp"
#include "statchem/molib/chain.hpp"
#include "statchem/molib/element.hpp"
#include "statchem/molib/grid.hpp"
#include "statchem/molib/it.hpp"
#include "statchem/molib/residue.hpp"

namespace statchem {

namespace molib {
class Chain;
class Residue;
class Atom;
class Assembly;

class Model : public template_map_container<Chain, Model, Assembly, char> {
    int __number;
    std::map<std::pair<int, Residue::res_tuple2>,
             std::map<Residue::res_tuple2, Residue::res_tuple2>>
        __remarks;
    Fragmenter::Fragment::Vec __rigid;

   public:
    Model(int number) : __number(number) {}
    Model(const Model& rhs) : __number(rhs.__number), __remarks(rhs.__remarks) {
        for (auto& chain : rhs) {
            dbgmsg("Copy constructor : model");
            add(new Chain(chain));
        }
    }
    Fragmenter::Fragment::Vec& get_rigid() { return __rigid; }
    const Fragmenter::Fragment::Vec& get_rigid() const { return __rigid; }
    void set_rigid(const Fragmenter::Fragment::Vec& rigid) { __rigid = rigid; }
    void init_bio(const Model& model_asym, const geometry::Matrix& matrix,
                  const std::set<char>& chains);
    void rotate(const geometry::Matrix& rota, const bool inverse = false);
    Chain& add(Chain* chain) {
        return this->aadd(chain->chain_id(), chain, this);
    }
    void set_remark(
        const int remark_number, const Residue::res_tuple2& ligand,
        std::pair<const Residue::res_tuple2&, const Residue::res_tuple2&>
            rpair) {
        this->__remarks[make_pair(remark_number, ligand)].insert(
            make_pair(rpair.first, rpair.second));
    }
    Model& regenerate_bonds(const Model&);
    void set_number(const int& i) { __number = i; }
    Chain& chain(const char chain_id) const { return this->element(chain_id); }
    bool has_chain(const char chain_id) const {
        return this->has_element(chain_id);
    }
    int number() const { return __number; }

    bool remarks(const int remark_number, const Residue::res_tuple2 ligand,
                 std::map<Residue::res_tuple2, Residue::res_tuple2>& r) {
        auto i = __remarks.find(make_pair(remark_number, ligand));
        if (i == __remarks.end()) return false;
        r = i->second;
        return true;
    }

    Atom::Vec get_atoms(
        const std::string& chain_ids = "",
        const Residue::res_type& rest = Residue::res_type::notassigned) const;
    Model& erase_properties() {
        for (auto& chain : *this) chain.erase_properties();
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& stream, const Model& m);
};

}  // molib
}

#endif
