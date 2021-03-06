/* This is score.hpp and is part of StatChemLIB
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

#ifndef SCORE_H
#define SCORE_H
#include <cmath>
#include <exception>
#include <map>
#include <set>
#include <typeinfo>
#include "statchem/geometry/geometry.hpp"
#include "statchem/helper/array1d.hpp"
#include "statchem/helper/error.hpp"
#include "statchem/molib/atom.hpp"

namespace statchem {

namespace molib {
class Molecule;
class Molecules;
}  // namespace molib

namespace score {

typedef std::pair<int, int> pair_of_ints;
typedef std::map<pair_of_ints, std::vector<double>> AtomPairValues;
typedef std::map<pair_of_ints, double> AtomPairSum;

struct AtomicDistributions {
    AtomPairValues values;
    double step_in_file;
    size_t max_distance;
    AtomicDistributions(const std::string& filename);
};

class Score {
   protected:
    AtomPairValues __gij_of_r_numerator;
    AtomPairValues __energies_scoring;  // scoring function

    AtomPairSum __sum_gij_of_r_numerator;

    std::vector<double> __gij_of_r_bin_range_sum, __bin_range_sum;
    double __total_quantity;
    std::set<pair_of_ints> __prot_lig_pairs;
    std::set<pair_of_ints> __avail_prot_lig;

    const double __eps;
    const std::string __ref_state, __comp, __rad_or_raw;
    const double __dist_cutoff;

    double __step_in_file;

    double __energy_mean(const pair_of_ints&, const double&);
    double __energy_cumulative(const pair_of_ints&, const double&);

    int __get_index(const double d) const {
        return (int)floor((d + 0.0000000001) / (double)__step_in_file);
    }
    double __get_lower_bound(const int idx) const {
        return (double)idx * (double)__step_in_file;
    }

   public:
    Score(const std::string& ref_state, const std::string& comp,
          const std::string& rad_or_raw, const double& dist_cutoff)
        : __total_quantity(0),
          __eps(0.0000001),
          __ref_state(ref_state),
          __comp(comp),
          __rad_or_raw(rad_or_raw),
          __dist_cutoff(dist_cutoff),
          __step_in_file(-1) {}

    double non_bonded_energy(const molib::Atom::Grid& gridrec,
                             const molib::Molecule&) const;
    double non_bonded_energy(const molib::Atom::Grid& gridrec,
                             const molib::Atom::Vec& atoms,
                             const geometry::Point::Vec& crds) const;

    Array1d<double> compute_energy(
        const molib::Atom::Grid& gridrec, const geometry::Coordinate& crd,
        const std::set<int>& ligand_atom_types) const;

    double get_dist_cutoff() const { return __dist_cutoff; }

    Score& define_composition(const std::set<int>& receptor_idatm_types,
                              const std::set<int>& ligand_idatm_types);

    Score& process_distributions(const AtomicDistributions& distributions);
    Score& process_distributions(const std::string& distributions_file);
    Score& compile_scoring_function();

    friend std::ostream& operator<<(std::ostream& stream,
                                    const std::vector<double>& energy);
    friend std::ostream& operator<<(std::ostream& stream,
                                    const AtomPairValues& energies);
    friend std::ostream& operator<<(std::ostream& stream, const Score& score);
};
}  // namespace score
}  // namespace statchem

#endif
