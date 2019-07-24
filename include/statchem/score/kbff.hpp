/* This is kbff.hpp and is part of StatChemLIB
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

#ifndef KBFF_H
#define KBFF_H
#include "statchem/score/score.hpp"

namespace statchem {

namespace score {
class KBFF : public Score {
    AtomPairValues __energies;  // objective function
    double __step_non_bond;
    std::set<pair_of_ints> __unavailible;

   public:
    KBFF(const std::string& ref_state, const std::string& comp,
         const std::string& rad_or_raw, const double& dist_cutoff,
         const double& step_non_bond)
        : Score(ref_state, comp, rad_or_raw, dist_cutoff),
          __step_non_bond(step_non_bond) {}

    double get_step_nonbond() const { return __step_non_bond; }
    const AtomPairValues& get_energies() const { return __energies; }

    const std::set<pair_of_ints>& get_unavailible() const { return __unavailible; }
    bool is_availible(int idatm1, int idatm2);

    KBFF& compile_objective_function();
    KBFF& parse_objective_function(const std::string& obj_dir,
                                   const double scale_non_bond,
                                   const size_t max_step);
    KBFF& output_objective_function(const std::string& obj_dir);

    friend std::ostream& operator<<(std::ostream& stream, const KBFF& kbff);
};
}
}

#endif
