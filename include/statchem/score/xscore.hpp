#ifndef XSCORE_H
#define XSCORE_H

#include <vector>
#include "statchem/molib/molecule.hpp"

namespace statchem {

namespace score {

std::array<double, 5> vina_xscore(const molib::Atom::Grid& gridrec,
                                  const molib::Atom::Vec& atoms);
}
}

#endif
