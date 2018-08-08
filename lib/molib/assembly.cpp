#include "statchem/molib/assembly.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "statchem/geometry/geometry.hpp"
#include "statchem/molib/atom.hpp"
#include "statchem/molib/chain.hpp"
#include "statchem/molib/model.hpp"
#include "statchem/molib/residue.hpp"
using namespace std;

namespace statchem {
namespace molib {
ostream& operator<<(ostream& stream, const Assembly& a) {
    stream << setw(6) << left << "REMARK   6 " << a.name() << " " << a.number()
           << endl;
    for (auto& model : a) {
        stream << model;
    }
    stream << "REMARK   6 END" << endl;
    return stream;
}

void Assembly::init_bio(const Assembly& asym,
                        map<int, geometry::Matrix>& matrices,
                        const set<char>& chains) {
    for (auto& i : matrices) {
        int matrix_number = i.first;
        geometry::Matrix& matrix = i.second;
        Model& last = add(new Model(matrix_number));
        last.init_bio(asym.first(), matrix,
                      chains);  // asymmetric unit has just one model - the 1-st
    }
}

void Assembly::rotate(const geometry::Matrix& rota, const bool inverse) {
    for (auto& model : *this) {
        model.rotate(rota, inverse);
    }
}

Atom::Vec Assembly::get_atoms(const string& chain_ids,
                              const Residue::res_type& rest,
                              const int model_number) const {
    Atom::Vec atoms;
    for (auto& model : *this) {
        if (model_number == -1 || model.number() == model_number) {
            auto ret = model.get_atoms(chain_ids, rest);
            atoms.insert(atoms.end(), ret.begin(), ret.end());
        }
    }
    return atoms;
}
}
}
