#ifndef BONDTYPE_H
#define BONDTYPE_H
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "statchem/molib/it.hpp"
#include "statchem/molib/atom.hpp"
#include "statchem/molib/bond.hpp"

namespace statchem {

namespace molib {
class Atom;
struct AtomParams {
    int val;
    int aps;
    int con;
};
typedef std::map<Atom*, AtomParams> ValenceState;
typedef std::vector<ValenceState> ValenceStateVec;
typedef std::map<Bond*, int> BondToOrder;

class BondOrderError : public Error {
   public:
    BondOrderError(const std::string& msg) : Error(msg) {}
};

class BondOrder {
    static ValenceStateVec __create_valence_states(
        const Atom::Vec& atoms, const int max_valence_states);
    static void __dfs(const int level, const int sum, const int tps,
                      const std::vector<std::vector<AtomParams>>& V,
                      std::vector<AtomParams>& Q,
                      std::vector<std::vector<AtomParams>>& valence_states,
                      const int max_valence_states);
    static bool __discrepancy(const ValenceState& valence_state);
    // Windows defines the macro __success !
    static bool __my_success(const ValenceState& valence_state);
    static bool __basic_rules(ValenceState& valence_state,
                              BondToOrder& bond_orders);
    static void __trial_error(ValenceState& valence_state,
                              BondToOrder& bond_orders);
    static Bond& __get_first_unassigned_bond(const ValenceState& valence_state,
                                             BondToOrder& bond_orders);

   public:
    static void compute_rotatable_bonds(const Atom::Vec& atoms);
    static void compute_bond_order(const Atom::Vec& atoms);
    static void compute_bond_gaff_type(const Atom::Vec& atoms);
};

void compute_chirality(const Atom::Vec& bonds);

std::ostream& operator<<(std::ostream& stream,
                         const ValenceState& valence_state);
std::ostream& operator<<(std::ostream& stream,
                         const ValenceStateVec& valence_states);
std::ostream& operator<<(std::ostream& os, const BondToOrder& bond_orders);
}
}

#endif
