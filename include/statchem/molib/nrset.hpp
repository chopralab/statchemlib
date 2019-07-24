/* This is nrset.hpp and is part of StatChemLIB
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

#ifndef NRSET_H
#define NRSET_H
#include "statchem/molib/it.hpp"
#include "statchem/molib/molecule.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/molib/residue.hpp"

#include <random>

namespace statchem {

namespace molib {

class NRset : public template_map_container<Molecules, NRset, NRset> {
   public:
    Molecules& add(Molecules* m) { return this->aadd(m, this); }
    NRset& erase_properties() {
        for (auto& molecules : *this) molecules.erase_properties();
        return *this;
    }
    Molecule::Vec get_molecules(const Residue::res_type& rest) const;
    Atom::Vec get_atoms(
        const std::string& chain_ids = "",
        const Residue::res_type& rest = Residue::res_type::notassigned,
        const int model_number = -1) const;
    void jiggle(std::mt19937 rng);
    friend std::ostream& operator<<(std::ostream& stream, const NRset& m);
};

}  // molib
}

#endif
