/* This is interpolation.hpp and is part of StatChemLIB
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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <vector>
#include <cstdlib>

namespace statchem {

namespace score {
namespace Interpolation {

std::vector<double> derivative(const std::vector<double>& y, const double step);
std::vector<double> interpolate(const std::vector<double>& dataX,
                                const std::vector<double>& dataY,
                                const double step);

struct BSplineFitWorkspace;

class BSplineFit {
    const size_t n;
    const size_t ncoeffs2;
    const size_t nbreak;

    BSplineFitWorkspace* bsplinewsp;

   public:
    BSplineFit(const std::vector<double>& dataX,
               const std::vector<double>& dataY, const size_t k = 4,
               const size_t ncoeff = 24);
    virtual ~BSplineFit();
    std::vector<double> interpolate_bspline(const double start,
                                            const double end,
                                            const double step);
};
}
}
}

#endif
