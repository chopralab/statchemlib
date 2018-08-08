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
