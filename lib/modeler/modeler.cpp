/* This is modeler.cpp and is part of StatChemLIB
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

#include "statchem/modeler/modeler.hpp"
#include <stdlib.h>
#include <time.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <regex>
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/benchmark.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/error.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/modeler/forcefield.hpp"
#include "statchem/modeler/systemtopology.hpp"
#include "statchem/modeler/topology.hpp"
#include "statchem/score/score.hpp"

#include <openmm/Units.h>
#include <openmm/Vec3.h>

using namespace std;

namespace statchem {
namespace OMMIface {
ostream& operator<<(ostream& os, const vector<OpenMM::Vec3>& positions) {
    os << "SIZE OF POSITIONS = " << positions.size() << endl;
    os << "ELEMENTS :" << endl;
    for (auto& v : positions)
        os << setprecision(8) << fixed << v[0] << " " << v[1] << " " << v[2]
           << endl;
    return os;
}

Modeler::Modeler(const ForceField& ffield, const string& fftype, double scale,
                 double tolerance, int max_iterations, bool use_constraints,
                 double step_size_in_fs, double temperature, double friction,
                 double cutoff)
    : __ffield(&ffield),
      __fftype(fftype),
      __scale(scale),
      __tolerance(tolerance),
      __max_iterations(max_iterations),
      __use_constraints(use_constraints),
      __step_size_in_ps(step_size_in_fs * OpenMM::PsPerFs),
      __temperature(temperature),
      __friction(friction),
      __cutoff(cutoff) {}

void Modeler::mask(const molib::Atom::Vec& atoms) {
    dbgmsg("Masking atoms " << atoms);
    __system_topology.mask(__topology, atoms);
}

void Modeler::unmask(const molib::Atom::Vec& atoms) {
    dbgmsg("Unmasking atoms " << atoms);
    __system_topology.unmask(__topology, atoms);
}

void Modeler::add_topology(const molib::Atom::Vec& atoms) {
    __topology.add_topology(atoms, *__ffield);
    __positions.resize(
        __topology.atoms.size());  // as many positions as there are atoms
}

void Modeler::add_crds(const molib::Atom::Vec& atoms,
                       const geometry::Point::Vec& crds) {
    for (size_t i = 0; i < atoms.size(); ++i) {
        int idx = __topology.get_index(*atoms[i]);
        __positions[idx] = crds[i];
    }
}

void Modeler::add_random_crds(const molib::Atom::Vec& atoms) {
    srand(time(NULL));
    for (size_t i = 0; i < atoms.size(); ++i) {
        int idx = __topology.get_index(*atoms[i]);
        const double x = rand() % 100 + 1;
        const double y = rand() % 100 + 1;
        const double z = rand() % 100 + 1;
        __positions[idx] = geometry::Point(x, y, z);
    }
}

/**
 * Changes coordinates of atoms
 */

geometry::Point::Vec Modeler::get_state(const molib::Atom::Vec& atoms) {
    const auto positions_in_nm = __system_topology.get_positions_in_nm();
    geometry::Point::Vec crds;
    crds.reserve(atoms.size());
    for (size_t i = 0; i < atoms.size(); ++i) {
        int idx = __topology.get_index(*atoms[i]);
        crds.emplace_back(
            geometry::Point(positions_in_nm[idx].x() * OpenMM::AngstromsPerNm,
                            positions_in_nm[idx].y() * OpenMM::AngstromsPerNm,
                            positions_in_nm[idx].z() * OpenMM::AngstromsPerNm));
    }

    return crds;
}

void Modeler::minimize_state() {
    Benchmark bench;
    __system_topology.minimize(__tolerance, __max_iterations);
    log_benchmark << "Minimization took " << bench.seconds_from_start()
                  << " wallclock seconds"
                  << "\n";
}

void Modeler::dynamics() {
    log_step << "Running " << __step_size_in_ps * __dynamics_steps
             << "ps simulation at " << __temperature << "K" << endl;
    __system_topology.dynamics(__dynamics_steps);
}

void Modeler::init_openmm_positions() {
#ifdef STATCHEM_DEBUG_MESSAGES
    dbgmsg("init_openmm_positions");
    for (auto& point : __positions) dbgmsg(point);
#endif
    __system_topology.init_positions(__positions);
}

void Modeler::init_openmm(const std::string& platform,
                          const std::string& precision,
                          const std::string& accelerators,
                          SystemTopology::integrator_type type) {
    __system_topology.set_forcefield(*__ffield);
    __system_topology.init_particles(__topology);
    __system_topology.init_bonded(__topology, __use_constraints);

    if (__fftype == "kb") {
        __system_topology.init_knowledge_based_force_3d(__topology, __scale,
                                                     __cutoff);
    } else if (__fftype == "phy") {
        __system_topology.init_physics_based_force(__topology);
    } else if (__fftype == "none") {
        // Do nothing
    } else {
        throw Error("die : unsupported forcefield");
    }

    __system_topology.init_integrator(type, __step_size_in_ps, __temperature,
                                      __friction);
    __system_topology.init_platform(platform, precision, accelerators);
}

double Modeler::potential_energy() {
    return __system_topology.get_potential_energy();
}
}  // namespace OMMIface
}  // namespace statchem
