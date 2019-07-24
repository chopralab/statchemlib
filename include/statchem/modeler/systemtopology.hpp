/* This is systemtopology.hpp and is part of StatChemLIB
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

#ifndef SYSTEMTOPOLOGY_H
#define SYSTEMTOPOLOGY_H
#include <map>
#include <set>
#include <string>
#include <vector>
#include "statchem/geometry/geometry.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/modeler/topology.hpp"
#include "statchem/molib/molecule.hpp"

namespace OpenMM {
class System;
class Integrator;
class Context;
class HarmonicAngleForce;
class HarmonicBondForce;
class PeriodicTorsionForce;
class CustomNonbondedForce;
}  // namespace OpenMM

namespace statchem {

namespace molib {
class Atom;
class Molecule;
}  // namespace molib

namespace OMMIface {
struct ForceField;

class SystemTopology {
   public:
    enum integrator_type {
        none,
        verlet,
        langevin,
        brownian,
    };

   private:
    #define num_checkpoints 5
    int checkpoint_num;
    
    OpenMM::System* system;
    OpenMM::Integrator* integrator;
    OpenMM::Context* context;
    OpenMM::CustomNonbondedForce* forcefield;

    integrator_type __integrator_used;
    int __thermostat_idx;

    OpenMM::HarmonicBondForce* bondStretch;
    OpenMM::HarmonicAngleForce* bondBend;
    OpenMM::PeriodicTorsionForce* bondTorsion;

    const ForceField* __ffield;

    std::vector<int> __kbforce_idx;
    std::vector<bool> masked;
    std::vector<double> masses;

    class AtomPoint {
       private:
        const geometry::Point __crd;
        molib::Atom& __atom;

       public:
        AtomPoint(const geometry::Point& crd, molib::Atom& atom)
            : __crd(crd), __atom(atom) {}
        const geometry::Point& crd() const { return __crd; }
        molib::Atom& get_atom() { return __atom; }
        void distance(double) const {}  // just dummy : needed by grid

        typedef std::vector<std::unique_ptr<AtomPoint>> UPVec;
        typedef std::vector<AtomPoint*> PVec;
        typedef statchem::molib::Grid<AtomPoint> Grid;
    };

    struct ForceData {
        int force_idx, idx1, idx2, idx3, idx4;
        double length, angle;
        int periodicity;
        double phase, k;
    };
    std::vector<std::vector<ForceData>> bondStretchData, bondBendData,
        bondTorsionData;

    void retype_amber_protein_atom_to_gaff(const molib::Atom& atom, int& type);

   public:
    SystemTopology()
        : system(nullptr),
          integrator(nullptr),
          context(nullptr),
          __integrator_used(integrator_type::none),
          __thermostat_idx(-1) {}
    ~SystemTopology();
    static void loadPlugins(const std::string& extra_dir = "");
    void mask(Topology& topology, const molib::Atom::Vec& atoms);
    void unmask(Topology& topology, const molib::Atom::Vec& atoms);

    void mask_forces(const int atom_idx, const std::set<int>& substruct);
    void unmask_forces(const int atom_idx, const std::set<int>& substruct);

    void init_integrator(SystemTopology::integrator_type type,
                         const double step_size_in_ps,
                         const double temperature_in_kelvin,
                         const double friction_in_per_ps);

    void init_platform(const std::string& platform,
                       const std::string& precision,
                       const std::string& accelerators);

    void init_particles(Topology& topology);
    void init_physics_based_force(Topology& topology);
    void init_knowledge_based_force(Topology& topology, double scale,
                                    double cutoff);
    void init_knowledge_based_force_3d(Topology& topology,
                                       double scale, double cutoff);
    void init_bonded(Topology& topology, const bool use_constraints);
    void init_positions(const geometry::Point::Vec& crds);

    void update_thermostat(const double temperature_in_kelvin,
                           const double collision_frequency);

    geometry::Point::Vec get_positions_in_nm();
    geometry::Point::Vec get_forces();

    double get_potential_energy();
    double get_kinetic_energy();

    void set_temperature();
    void set_box_vector();

    void load_checkpoint(const std::string& checkpoint);
    void save_checkpoint();
    void save_checkpoint_candock(int x, int y);

    // Print kinetic, potential, and total energies to stderr
    void print_energies();

    void print_box_vector_size();

    void minimize(const double tolerance, const int max_iterations);
    void dynamics(const int steps);
    double get_energies();
    void set_forcefield(const ForceField& ffield) { __ffield = &ffield; }
};
}  // namespace OMMIface
}  // namespace statchem

#endif
