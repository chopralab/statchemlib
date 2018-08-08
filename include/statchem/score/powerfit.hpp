#ifndef POWERFIT_H
#define POWERFIT_H

#include <tuple>
#include "statchem/helper/error.hpp"

namespace statchem {

namespace score {

class InterpolationError : public Error {
   public:
    InterpolationError(const std::string& msg) : Error(msg) {}
};

std::tuple<double, double, double> fit_power_function(
    std::vector<double> x, std::vector<double> y, size_t power,
    std::pair<double, double> guess, size_t max_iter);

std::tuple<double, double, double> fit_range_power_function(
    std::vector<double> x, std::vector<double> y);
std::tuple<double, double, double> fit_range_power_function_fast(
    std::vector<double> x, std::vector<double> y);
}
}

#endif
