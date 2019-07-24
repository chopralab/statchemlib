/* This is powerfit.hpp and is part of StatChemLIB
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
