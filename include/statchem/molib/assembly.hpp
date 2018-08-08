#ifndef ASSEMBLY_H
#define ASSEMBLY_H
#include "statchem/geometry/geometry.hpp"
#include "statchem/geometry/matrix.hpp"
#include "statchem/molib/atom.hpp"
#include "statchem/molib/chain.hpp"
#include "statchem/molib/element.hpp"
#include "statchem/molib/grid.hpp"
#include "statchem/molib/it.hpp"
#include "statchem/molib/model.hpp"
#include "statchem/molib/residue.hpp"

namespace statchem {

namespace molib {
class Chain;
class Residue;
class Atom;
class Molecule;

class Assembly : public template_map_container<Model, Assembly, Molecule> {
    int __number;
    std::string __name;  // ASYMMETRIC UNIT OR BIOLOGICAL ASSEMBLY
   public:
    Assembly(int number, const std::string name = "ASYMMETRIC UNIT")
        : __number(number), __name(name) {}
    Assembly(const Assembly& rhs) : __number(rhs.__number), __name(rhs.__name) {
        for (auto& model : rhs) {
            dbgmsg("Copy constructor : assembly");
            add(new Model(model));
        }
    }

    void init_bio(const Assembly& asym,
                  std::map<int, geometry::Matrix>& matrices,
                  const std::set<char>& chains);

    void set_number(int number) { __number = number; }
    void set_name(const std::string& name) { __name = name; }
    int number() const { return __number; }
    std::string name() const { return __name; }
    Model& add(Model* m) { return this->aadd(m->number(), m, this); }

    void rotate(const geometry::Matrix& rota, const bool inverse = false);

    Atom::Vec get_atoms(
        const std::string& chain_ids = "",
        const Residue::res_type& rest = Residue::res_type::notassigned,
        const int model_number = -1) const;
    Assembly& erase_properties() {
        for (auto& model : *this) model.erase_properties();
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Assembly& a);
};

}  // molib
}

#endif
