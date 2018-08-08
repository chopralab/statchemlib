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
