#include "statchem/geometry/coordinate.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "statchem/geometry/matrix.hpp"
using namespace std;

namespace statchem {
namespace geometry {

unique_ptr<Coordinate> Coordinate::rotate(const Matrix& matrix) {
    unique_ptr<Coordinate> c(new Coordinate(*this));
    matrix.rotate(*c);
    return c;
}

unique_ptr<Coordinate> Coordinate::inverse_rotate(const Matrix& matrix) {
    unique_ptr<Coordinate> c(new Coordinate(*this));
    matrix.inverse_rotate(*c);
    return c;
}
}
}
