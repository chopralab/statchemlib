/* This is kabsch.hpp and is part of StatChemLIB
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

#ifndef KABSCH_H
#define KABSCH_H
#include "statchem/geometry/coordinate.hpp"
#include "statchem/geometry/matrix.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/error.hpp"
#include <memory>

namespace statchem {

class Kabsch {

    struct KabschPrivate;
    std::unique_ptr<KabschPrivate> __private;

    int __counter, __sz;

   public:
    Kabsch(const int sz = 0);
    void resize(const int sz);
    void clear();
    void add_vertex(const geometry::Coordinate& c,
                    const geometry::Coordinate& d);
    void superimpose();
    //geometry::Matrix get_rota() const { return geometry::Matrix(__U, __t); }
};
}

#endif  // KABSCH_H
