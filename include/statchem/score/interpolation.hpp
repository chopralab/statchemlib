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
