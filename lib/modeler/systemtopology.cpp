/* This is systemtopology.cpp and is part of StatChemLIB
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

#include "statchem/modeler/systemtopology.hpp"
#include "statchem/modeler/modeler.hpp"

#include <openmm/AndersenThermostat.h>
#include <openmm/BrownianIntegrator.h>
#include <openmm/CustomNonbondedForce.h>
#include <openmm/HarmonicAngleForce.h>
#include <openmm/HarmonicBondForce.h>
#include <openmm/LangevinIntegrator.h>
#include <openmm/LocalEnergyMinimizer.h>
#include <openmm/NonbondedForce.h>
#include <openmm/PeriodicTorsionForce.h>
#include <openmm/Platform.h>
#include <openmm/Units.h>
#include <openmm/Vec3.h>
#include <openmm/VerletIntegrator.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <regex>
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/error.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/modeler/forcefield.hpp"
#include "statchem/modeler/topology.hpp"

#include <chrono>
#include <map>
#include <utility>  // std::pair, std::make_pair
#include <iostream>
#include <fstream>


using namespace std;
using namespace std::chrono;

namespace statchem {
namespace OMMIface {

SystemTopology::~SystemTopology() {
    dbgmsg("calling destructor of SystemTopology");
    delete context;
    // delete integrator;
    delete system;
}

void SystemTopology::loadPlugins(const std::string& extra_dir) {
    try {
        // Load all available OpenMM plugins from their default location.

        dbgmsg("before loading plugins");
        OpenMM::Platform::loadPluginsFromDirectory(
            OpenMM::Platform::getDefaultPluginsDirectory());

        if (!extra_dir.empty()) {
            OpenMM::Platform::loadPluginsFromDirectory(extra_dir);
        }

        dbgmsg("after loading plugins");
    } catch (const std::exception& e) {
        log_error << "The crash is " << e.what() << endl;
    }
}

void SystemTopology::mask(Topology& topology, const molib::Atom::Vec& atoms) {
    set<int> substruct;

    for (auto& patom : atoms) {
        int idx = topology.get_index(*patom);
        dbgmsg("masking particle idx = " << idx);
        substruct.insert(idx);
    }

    for (auto& idx : substruct) {
        system->setParticleMass(idx, 0);
        masked[idx] = true;

        mask_forces(idx, substruct);
    }

    bondStretch->updateParametersInContext(*context);
    bondBend->updateParametersInContext(*context);
    bondTorsion->updateParametersInContext(*context);
}

void SystemTopology::unmask(Topology& topology, const molib::Atom::Vec& atoms) {
    set<int> substruct;

    for (auto& patom : atoms) {
        int idx = topology.get_index(*patom);
        substruct.insert(idx);
    }

    for (auto& idx : substruct) {
        dbgmsg("unmasking particle idx = " << idx << " mass = " << masses[idx]);
        system->setParticleMass(idx, masses[idx]);
        masked[idx] = false;

        unmask_forces(idx, substruct);
    }

    bondStretch->updateParametersInContext(*context);
    bondBend->updateParametersInContext(*context);
    bondTorsion->updateParametersInContext(*context);
}

void SystemTopology::mask_forces(const int atom_idx,
                                 const set<int>& substruct) {
    for (auto& data : bondStretchData[atom_idx]) {  // get all forces involving
                                                    // this atom's idx
        if (substruct.count(data.idx1) && substruct.count(data.idx2))
            bondStretch->setBondParameters(data.force_idx, data.idx1, data.idx2,
                                           data.length, 0.0);
    }

    for (auto& data :
         bondBendData[atom_idx]) {  // get all forces involving this atom's idx
        if (substruct.count(data.idx1) && substruct.count(data.idx2) &&
            substruct.count(data.idx3))
            bondBend->setAngleParameters(data.force_idx, data.idx1, data.idx2,
                                         data.idx3, data.angle, 0.0);
    }

    for (auto& data : bondTorsionData[atom_idx]) {  // get all forces involving
                                                    // this atom's idx
        if (substruct.count(data.idx1) && substruct.count(data.idx2) &&
            substruct.count(data.idx3) && substruct.count(data.idx4)) {
            int idx1, idx2, idx3, idx4;
            int periodicity;
            double phase, k;

            bondTorsion->getTorsionParameters(data.force_idx, idx1, idx2, idx3,
                                              idx4, periodicity, phase, k);
            dbgmsg("particles of force " << data.force_idx);
            dbgmsg(idx1 << " " << data.idx1);
            dbgmsg(idx2 << " " << data.idx2);
            dbgmsg(idx3 << " " << data.idx3);
            dbgmsg(idx4 << " " << data.idx4);
            dbgmsg(periodicity << " " << data.periodicity);

            if (idx1 != data.idx1) throw Error("die : particles changed");

            if (idx2 != data.idx2) throw Error("die : particles changed");

            if (idx3 != data.idx3) throw Error("die : particles changed");

            if (idx4 != data.idx4) throw Error("die : particles changed");

            bondTorsion->setTorsionParameters(
                data.force_idx, data.idx1, data.idx2, data.idx3, data.idx4,
                data.periodicity, data.phase, 0.0);
        }
    }
}

void SystemTopology::unmask_forces(const int atom_idx,
                                   const set<int>& substruct) {
    for (auto& data : bondStretchData[atom_idx]) {  // get all forces involving
                                                    // this atom's idx
        if (substruct.count(data.idx1) && substruct.count(data.idx2))
            bondStretch->setBondParameters(data.force_idx, data.idx1, data.idx2,
                                           data.length, data.k);
    }

    for (auto& data :
         bondBendData[atom_idx]) {  // get all forces involving this atom's idx
        if (substruct.count(data.idx1) && substruct.count(data.idx2) &&
            substruct.count(data.idx3))
            bondBend->setAngleParameters(data.force_idx, data.idx1, data.idx2,
                                         data.idx3, data.angle, data.k);
    }

    for (auto& data : bondTorsionData[atom_idx]) {  // get all forces involving
                                                    // this atom's idx
        if (substruct.count(data.idx1) && substruct.count(data.idx2) &&
            substruct.count(data.idx3) && substruct.count(data.idx4)) {
            int idx1, idx2, idx3, idx4;
            int periodicity;
            double phase, k;

            bondTorsion->getTorsionParameters(data.force_idx, idx1, idx2, idx3,
                                              idx4, periodicity, phase, k);
            dbgmsg("particles of force " << data.force_idx);
            dbgmsg(idx1 << " " << data.idx1);
            dbgmsg(idx2 << " " << data.idx2);
            dbgmsg(idx3 << " " << data.idx3);
            dbgmsg(idx4 << " " << data.idx4);
            dbgmsg(periodicity << " " << data.periodicity);

            if (idx1 != data.idx1) throw Error("die : particles changed");

            if (idx2 != data.idx2) throw Error("die : particles changed");

            if (idx3 != data.idx3) throw Error("die : particles changed");

            if (idx4 != data.idx4) throw Error("die : particles changed");

            bondTorsion->setTorsionParameters(
                data.force_idx, data.idx1, data.idx2, data.idx3, data.idx4,
                data.periodicity, data.phase, data.k);
        }
    }
}

// Choose an Integrator for advancing time, and a Context connecting the
// System with the Integrator for simulation. Let the Context choose the
// best available Platform. Initialize the configuration from the default
// positions we collected above. Initial velocities will be zero.

void SystemTopology::init_integrator(SystemTopology::integrator_type type,
                                     const double step_size_in_ps,
                                     const double temperature_in_kelvin,
                                     const double friction_in_per_ps) {
    if (integrator != nullptr) throw Error("Integrator already initialized");

    __integrator_used = type;

    switch (__integrator_used) {
        case verlet:
            __thermostat_idx = system->addForce(new OpenMM::AndersenThermostat(
                temperature_in_kelvin, friction_in_per_ps));
        // Fall through
        case none:
        default:
            integrator = new OpenMM::VerletIntegrator(step_size_in_ps);
            break;
        case langevin:
            integrator = new OpenMM::LangevinIntegrator(
                temperature_in_kelvin, friction_in_per_ps, step_size_in_ps);
            break;
        case brownian:
            integrator = new OpenMM::BrownianIntegrator(
                temperature_in_kelvin, friction_in_per_ps, step_size_in_ps);
            break;
    }
}

void SystemTopology::init_platform(const std::string& platform,
                                   const std::string& precision,
                                   const std::string& accelerators) {
    map<string, string> properties;

    if (platform == "CUDA" || platform == "OpenCL") {
        properties["Precision"] = precision;
        properties["DeviceIndex"] = accelerators;
    }

    // Available Platforms: Reference, CPU, CUDA, and OpenCL
    context = new OpenMM::Context(*system, *integrator,
                                  OpenMM::Platform::getPlatformByName(platform),
                                  properties);

    dbgmsg("Using OpenMM platform " << context->getPlatform().getName());
}

void SystemTopology::init_particles(Topology& topology) {
    int warn = 0;

    system = new OpenMM::System();

    dbgmsg("topology.atoms.size() = " << topology.atoms.size());

    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD
    //  units!).
    //  (3) Collect default positions for initializing the simulation later.

    for (auto& patom : topology.atoms) {
        const molib::Atom& atom = *patom;
        const int type = topology.get_type(atom);

        try {
            const ForceField::AtomType& atype = __ffield->get_atom_type(type);
            system->addParticle(atype.mass);

            masses.push_back(atype.mass);
            masked.push_back(false);  // at the start unmask everything
            dbgmsg("add particle type = "
                   << type << " crd = " << atom.crd()
                   << " mass = " << atype.mass << " charge = " << atype.charge
                   << " sigma = " << atype.sigma << " epsilon = "
                   << atype.epsilon << " representing atom = " << atom
                   << " at index = " << topology.get_index(atom));
        } catch (ParameterError& e) {
            log_error << e.what() << " (" << ++warn << ")" << endl;
        }
    }

    if (warn > 0) {
        throw Error("die : missing parameters detected");
    }
}

void SystemTopology::update_thermostat(const double temperature_in_kelvin,
                                       const double collision_frequency) {
    switch (__integrator_used) {
        case verlet:
            dynamic_cast<OpenMM::AndersenThermostat&>(
                system->getForce(__thermostat_idx))
                .setDefaultTemperature(temperature_in_kelvin);
            dynamic_cast<OpenMM::AndersenThermostat&>(
                system->getForce(__thermostat_idx))
                .setDefaultCollisionFrequency(collision_frequency);
            break;
        case langevin:
            dynamic_cast<OpenMM::LangevinIntegrator*>(integrator)
                ->setTemperature(temperature_in_kelvin);
            dynamic_cast<OpenMM::LangevinIntegrator*>(integrator)
                ->setFriction(collision_frequency);
            break;
        case brownian:
            dynamic_cast<OpenMM::BrownianIntegrator*>(integrator)
                ->setTemperature(temperature_in_kelvin);
            dynamic_cast<OpenMM::BrownianIntegrator*>(integrator)
                ->setFriction(collision_frequency);
            break;
        case none:
        default:
            throw Error("No integrator is being used!");
            break;
    }
}

void SystemTopology::init_physics_based_force(Topology& topology) {
    int warn = 0;

    OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
    nonbond->setNonbondedMethod(OpenMM::NonbondedForce::NonbondedMethod::PME);
    nonbond->setCutoffDistance(2.99);
    system->addForce(nonbond);
    system->setDefaultPeriodicBoxVectors(
        OpenMM::Vec3(6, 0, 0), OpenMM::Vec3(0, 6, 0), OpenMM::Vec3(0, 0, 6));

    for (auto& patom : topology.atoms) {
        const molib::Atom& atom = *patom;
        const int type = topology.get_type(atom);

        try {
            const ForceField::AtomType& atype = __ffield->get_atom_type(type);

            nonbond->addParticle(atype.charge, atype.sigma, atype.epsilon);

            dbgmsg("add particle type = "
                   << type << " crd = " << atom.crd()
                   << " mass = " << atype.mass << " charge = " << atype.charge
                   << " sigma = " << atype.sigma << " epsilon = "
                   << atype.epsilon << " representing atom = " << atom
                   << " at index = " << topology.get_index(atom));
        } catch (ParameterError& e) {
            log_error << e.what() << " (" << ++warn << ")" << endl;
        }
    }

    vector<pair<int, int>> bondPairs;

    for (auto& bond : topology.bonds) {
        dbgmsg("checkpoint0");
        const molib::Atom& atom1 = *bond.first;
        dbgmsg(atom1);
        const molib::Atom& atom2 = *bond.second;
        dbgmsg(atom2);
        const int idx1 = topology.get_index(atom1);
        dbgmsg(idx1);
        const int idx2 = topology.get_index(atom2);
        dbgmsg(idx2);
        bondPairs.push_back({idx1, idx2});
    }

    dbgmsg("checkpoint3");
    // Exclude 1-2, 1-3 bonded atoms from nonbonded forces, and scale down 1-4
    // bonded atoms.
    nonbond->createExceptionsFromBonds(bondPairs, __ffield->coulomb14scale,
                                       __ffield->lj14scale);

    if (warn > 0) {
        throw Error("die : missing parameters detected");
    }
}

void SystemTopology::init_knowledge_based_force(Topology& topology,
                                                double scale, double cutoff) {
    while (__kbforce_idx.size() > 0) {
        system->removeForce(__kbforce_idx.back());
        __kbforce_idx.pop_back();
    }

    std::set<int> used_atom_types;

    int position = 0;

    // First: atom type
    // Second: position in OpenMM
    multimap<int, int> idatm_mapping;

    // Map the idatm_type to the position in OpenMM
    for (const auto& atom : topology.atoms) {
        used_atom_types.insert(atom->idatm_type());
        idatm_mapping.insert(pair<int, int>(atom->idatm_type(), position++));
    }

    try {
        vector<pair<int, int>> bondPairs;

        for (auto& bond : topology.bonds) {
            const molib::Atom& atom1 = *bond.first;
            const molib::Atom& atom2 = *bond.second;
            const int idx1 = topology.get_index(atom1);
            const int idx2 = topology.get_index(atom2);
            bondPairs.push_back({idx1, idx2});
        }

        for (auto idatm1 = used_atom_types.begin();
             idatm1 != used_atom_types.end(); idatm1++) {
            for (auto idatm2 = idatm1; idatm2 != used_atom_types.end();
                 idatm2++) {
                // Create new CustomNonbondedForce
                auto forcefield =
                    new OpenMM::CustomNonbondedForce("scale * kbpot(r)");

                forcefield->setNonbondedMethod(
                    OpenMM::CustomNonbondedForce::CutoffNonPeriodic);

                // Cutoff distance needs to be less than half the size of the
                // Periodic Box Size
                forcefield->setCutoffDistance(__ffield->kb_cutoff);
                forcefield->addGlobalParameter("scale", scale);

                auto& pot = __ffield->kb_force_type.at(*idatm1)
                                                   .at(*idatm2)
                                      .potential;

                forcefield->addTabulatedFunction(
                    "kbpot", new OpenMM::Continuous1DFunction(
                                 pot, 0, (pot.size() - 1) * __ffield->step));

                set<int> one, two;
                auto type_1_iter = idatm_mapping.find(*idatm1);
                auto type_2_iter = idatm_mapping.find(*idatm2);

                for (; type_1_iter != idatm_mapping.end(); ++type_1_iter)
                    one.insert(type_1_iter->second);

                for (; type_2_iter != idatm_mapping.end(); ++type_2_iter)
                    two.insert(type_2_iter->second);

                forcefield->addInteractionGroup(one, two);

                // Create a bunch of empty particle parameters
                for (const auto& atom : topology.atoms)
                    forcefield->addParticle(vector<double>());

                forcefield->createExclusionsFromBonds(bondPairs, 4);

                __kbforce_idx.push_back(system->addForce(forcefield));
                break;
            }
        }
    } catch (ParameterError& e) {
        cerr << e.what() << endl;
        cerr << "Exiting" << endl;
        exit(0);
    }
}

void SystemTopology::init_knowledge_based_force_3d(Topology& topology,
                                                double scale, double cutoff) {

    forcefield = new OpenMM::CustomNonbondedForce(
        "scale * kbpot( r, idatm1 , idatm2)");
    forcefield->setNonbondedMethod(
        OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    forcefield->setCutoffDistance(__ffield->kb_cutoff);
    forcefield->addGlobalParameter("scale", scale);
    forcefield->addPerParticleParameter("idatm");

    std::map<int, int> __idatm_to_internal;
    std::map<int, int> __internal_to_idatm;
    int num_types = 0;
    for (const auto& atom : topology.atoms) {
        if (!__idatm_to_internal.count(atom->idatm_type())) {
            __idatm_to_internal[atom->idatm_type()] = num_types;
            __internal_to_idatm[num_types] = atom->idatm_type();
            num_types++;
        }

        forcefield->addParticle(
            {static_cast<double>(__idatm_to_internal[atom->idatm_type()])});
    }

    vector<double> table;
    size_t xsize = 0, ysize = 0;

    for (size_t i = 0; i < __idatm_to_internal.size(); i++) {
        for (size_t j = 0; j < __idatm_to_internal.size(); j++) {
            try {
                const ForceField::KBType kb = __ffield->get_kb_force_type(
                    __internal_to_idatm[i], __internal_to_idatm[j]);

                if (xsize == 0) xsize = kb.potential.size();

                if (xsize != kb.potential.size())
                    throw Error("mismatching potential size: Got " +
                                std::to_string(kb.potential.size()) +
                                " wanted " + std::to_string(xsize)  +
                                " for " + std::string(help::idatm_unmask[__internal_to_idatm[i]]) +
                                " and " + std::string(help::idatm_unmask[__internal_to_idatm[j]]));

                ysize++;

                for (size_t i = 0; i < kb.potential.size(); i++)
                    table.push_back(kb.potential[i]);
            } catch (ParameterError& e) {
                cerr << e.what() << endl;
                cerr << "This is normal for atom types only present in the "
                        "ligand. Adding " << xsize << " zeros"
                     << endl;
                for (size_t i = 0; i < xsize; i++)
                    table.push_back(0.0);
            }
        }
    }

    try {
        forcefield->addTabulatedFunction(
            "kbpot",
            new OpenMM::Continuous3DFunction(xsize,
                                             __idatm_to_internal.size(),
                                             __idatm_to_internal.size(), table,
                                             0.0, (xsize - 1) * __ffield->step,
                                             0.0, __idatm_to_internal.size() - 1.0,
                                             0.0, __idatm_to_internal.size() - 1.0
                                            ));

        vector<pair<int, int>> bondPairs;

        for (auto& bond : topology.bonds) {
            const molib::Atom& atom1 = *bond.first;
            const molib::Atom& atom2 = *bond.second;
            const int idx1 = topology.get_index(atom1);
            const int idx2 = topology.get_index(atom2);
            bondPairs.push_back({idx1, idx2});
        }

        forcefield->createExclusionsFromBonds(bondPairs, 4);
    } catch (ParameterError& e) {
        cerr << e.what() << endl;
        cerr << "Exiting" << endl;
        exit(0);
    }
    system->addForce(forcefield);
}

void SystemTopology::retype_amber_protein_atom_to_gaff(const molib::Atom& atom,
                                                       int& type) {
    // Only retype protein atoms to gaff atoms
    if (atom.br().rest() != molib::Residue::protein) {
        // Its the N in the peptide bond needs to be N, not n3....
        if (atom.element() == molib::Element::N) {
            type = __ffield->gaff_name_to_type.at("n");
        }
    }

    if (atom.element() == molib::Element::C) {
        // C2 == 18 aka the carbon is in the carbonyl
        if (atom.idatm_type() == 18) type = __ffield->gaff_name_to_type.at("c");
        // Otherwise it's the alpha carbon
        else
            type = __ffield->gaff_name_to_type.at("c3");
    }

    // Its the N in the peptide bond
    if (atom.element() == molib::Element::N) {
        type = __ffield->gaff_name_to_type.at("n");
    }

    // Its the O in the peptide bond
    if (atom.element() == molib::Element::O) {
        type = __ffield->gaff_name_to_type.at("o");
    }
}

void SystemTopology::init_bonded(Topology& topology,
                                 const bool use_constraints) {
    int warn = 0;

    bondStretchData.resize(topology.atoms.size());
    bondBendData.resize(topology.atoms.size());
    bondTorsionData.resize(topology.atoms.size());

    bondStretch = new OpenMM::HarmonicBondForce();
    bondBend = new OpenMM::HarmonicAngleForce();
    bondTorsion = new OpenMM::PeriodicTorsionForce();

    system->addForce(bondStretch);
    system->addForce(bondBend);
    system->addForce(bondTorsion);

    dbgmsg("initializing openmm");

    // Process the bonds:
    //  (1) If we're using constraints, tell System about constrainable bonds;
    //      otherwise, tell HarmonicBondForce the bond stretch parameters
    //      (tricky units!).
    //  (2) Create a list of bonds for generating nonbond exclusions.

    int force_idx = 0;

    for (auto& bond : topology.bonds) {
        dbgmsg("checkpoint0");

        const molib::Atom& atom1 = *bond.first;
        dbgmsg(atom1);
        const molib::Atom& atom2 = *bond.second;
        dbgmsg(atom2);
        const int idx1 = topology.get_index(atom1);
        dbgmsg(idx1);
        const int idx2 = topology.get_index(atom2);
        dbgmsg(idx2);
        int type1 = topology.get_type(atom1);
        dbgmsg(type1);
        int type2 = topology.get_type(atom2);
        dbgmsg(type2);

        // Check for modified residues
        int number_protein = (atom1.br().rest() == molib::Residue::protein) +
                             (atom2.br().rest() == molib::Residue::protein);

        int number_not_set =
            (atom1.gaff_type() != "???") + (atom2.gaff_type() != "???");

        // It should both protein or non-protein, otherwise we fix it
        if (number_protein == 1 && number_not_set == 1) {
            retype_amber_protein_atom_to_gaff(atom1, type1);
            retype_amber_protein_atom_to_gaff(atom2, type2);
        }

        ForceField::BondType btype;

        try {
            btype = __ffield->get_bond_type(type1, type2);
        } catch (ParameterError& e) {
            log_warning << e.what()
                        << " (WARNINGS ARE NOT INCREASED) (using "
                           "default parameters for this bond)"
                        << endl;
            // if everything else fails just constrain at something reasonable
            btype = ForceField::BondType{atom1.get_bond(atom2).length(), 250000,
                                         false};
        }

        if (use_constraints &&
            btype.can_constrain) {  // Should we constrain C-H bonds?
            system->addConstraint(idx1, idx2, btype.length);
        } else {
            // Note factor of 2 for stiffness below because Amber specifies the
            // constant
            // as it is used in the harmonic energy term kx^2 with force 2kx;
            // OpenMM wants
            // it as used in the force term kx, with energy kx^2/2.
            bondStretch->addBond(idx1, idx2, btype.length, btype.k);
            dbgmsg("force_idx = " << force_idx << " idx1 = " << idx1
                                  << " idx2 = " << idx2 << " bond length = "
                                  << btype.length << " k = " << btype.k);
            bondStretchData[idx1].push_back(ForceData{
                force_idx, idx1, idx2, 0, 0, btype.length, 0, 0, 0, btype.k});
            bondStretchData[idx2].push_back(ForceData{
                force_idx, idx1, idx2, 0, 0, btype.length, 0, 0, 0, btype.k});
            ++force_idx;
        }
    }

    dbgmsg("checkpoint3");

    force_idx = 0;

    // Create the 1-2-3 bond angle harmonic terms.
    for (auto& angle : topology.angles) {
        dbgmsg("checkpoint4");
        const molib::Atom& atom1 = *get<0>(angle);
        const molib::Atom& atom2 = *get<1>(angle);
        const molib::Atom& atom3 = *get<2>(angle);
        const int idx1 = topology.get_index(atom1);
        const int idx2 = topology.get_index(atom2);
        dbgmsg("checkpoint5");
        const int idx3 = topology.get_index(atom3);
        int type1 = topology.get_type(atom1);
        int type2 = topology.get_type(atom2);
        int type3 = topology.get_type(atom3);
        dbgmsg("checkpoint6");

        // Check for modified residues
        int number_protein = (atom1.br().rest() == molib::Residue::protein) +
                             (atom2.br().rest() == molib::Residue::protein) +
                             (atom3.br().rest() == molib::Residue::protein);

        int number_not_set = (atom1.gaff_type() != "???") +
                             (atom2.gaff_type() != "???") +
                             (atom3.gaff_type() != "???");

        // It should all be protein or non-protein, otherwise we fix it
        if ((number_protein == 1 || number_protein == 2) &&
            (number_not_set == 1 || number_not_set == 2)) {
            retype_amber_protein_atom_to_gaff(atom1, type1);
            retype_amber_protein_atom_to_gaff(atom2, type2);
            retype_amber_protein_atom_to_gaff(atom3, type3);
        }

        ForceField::AngleType atype;

        try {
            dbgmsg("determining angle type between atoms : " << endl
                                                             << atom1 << endl
                                                             << atom2 << endl
                                                             << atom3);
            atype = __ffield->get_angle_type(type1, type2, type3);
        } catch (ParameterError& e) {
            log_warning << e.what()
                        << " (WARNINGS ARE NOT INCREASED) (using "
                           "default parameters for this angle)"
                        << endl;
            // if everything else fails just constrain at something reasonable
            atype = ForceField::AngleType{
                geometry::angle(atom1.crd(), atom2.crd(), atom3.crd()), 500};
        }

        dbgmsg("force_idx = " << force_idx << " idx1 = " << idx1
                              << " idx2 = " << idx2 << " idx3 = " << idx3
                              << " angle = " << atype.angle
                              << " k = " << atype.k);
        bondBend->addAngle(idx1, idx2, idx3, atype.angle, atype.k);
        bondBendData[idx1].push_back(ForceData{force_idx, idx1, idx2, idx3, 0,
                                               0, atype.angle, 0, 0, atype.k});
        bondBendData[idx2].push_back(ForceData{force_idx, idx1, idx2, idx3, 0,
                                               0, atype.angle, 0, 0, atype.k});
        bondBendData[idx3].push_back(ForceData{force_idx, idx1, idx2, idx3, 0,
                                               0, atype.angle, 0, 0, atype.k});

        ++force_idx;

        dbgmsg("checkpoint8");
    }

    force_idx = 0;

    // Create the 1-2-3-4 bond torsion (dihedral) terms.
    for (auto& dihedral : topology.dihedrals) {
        dbgmsg("checkpoint9");
        const molib::Atom& atom1 = *get<0>(dihedral);
        const molib::Atom& atom2 = *get<1>(dihedral);
        const molib::Atom& atom3 = *get<2>(dihedral);
        const molib::Atom& atom4 = *get<3>(dihedral);
        const int idx1 = topology.get_index(atom1);
        const int idx2 = topology.get_index(atom2);
        const int idx3 = topology.get_index(atom3);
        const int idx4 = topology.get_index(atom4);
        int type1 = topology.get_type(atom1);
        int type2 = topology.get_type(atom2);
        int type3 = topology.get_type(atom3);
        int type4 = topology.get_type(atom4);

        // Check for modified residues
        int number_protein = (atom1.br().rest() == molib::Residue::protein) +
                             (atom2.br().rest() == molib::Residue::protein) +
                             (atom3.br().rest() == molib::Residue::protein) +
                             (atom4.br().rest() == molib::Residue::protein);

        int number_not_set =
            (atom1.gaff_type() != "???") + (atom2.gaff_type() != "???") +
            (atom3.gaff_type() != "???") + (atom4.gaff_type() != "???");

        // It should all be protein or non-protein, otherwise we fix it
        if ((number_protein == 1 || number_protein == 2 ||
             number_protein == 3) &&
            (number_not_set == 1 || number_not_set == 2 ||
             number_not_set == 3)) {
            retype_amber_protein_atom_to_gaff(atom1, type1);
            retype_amber_protein_atom_to_gaff(atom2, type2);
            retype_amber_protein_atom_to_gaff(atom3, type3);
            retype_amber_protein_atom_to_gaff(atom4, type4);
        }

        try {
            const ForceField::TorsionTypeVec& v_ttype =
                __ffield->get_dihedral_type(type1, type2, type3,
                                            type4);  // cannot make it const ??
            for (auto& ttype : v_ttype) {
                dbgmsg("force_idx = "
                       << force_idx << " idx1 = " << idx1 << " idx2 = " << idx2
                       << " idx3 = " << idx3 << " idx4 = " << idx4
                       << " periodicity = " << ttype.periodicity
                       << " phase = " << ttype.phase << " k = " << ttype.k);
                bondTorsion->addTorsion(idx1, idx2, idx3, idx4,
                                        ttype.periodicity, ttype.phase,
                                        ttype.k);
                bondTorsionData[idx1].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});
                bondTorsionData[idx2].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});
                bondTorsionData[idx3].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});
                bondTorsionData[idx4].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});

                ++force_idx;
            }
        } catch (ParameterError& e) {
            log_warning << e.what()
                        << " (WARNINGS ARE NOT INCREASED) (using "
                           "default parameters for this dihedral)"
                        << endl;
            // if everything else fails just constrain at something reasonable
            // cout <<  geometry::dihedral(atom1.crd(), atom2.crd(),
            // atom3.crd(), atom4.crd()) << endl;
            // ForceField::TorsionType vtype = ForceField::TorsionType {1,
            // geometry::dihedral(atom1.crd(), atom2.crd(), atom3.crd(),
            // atom4.crd()), 0};
            // bondTorsion->addTorsion (idx1, idx2, idx3, idx4,
            // vtype.periodicity, vtype.phase, vtype.k);
            // bondTorsionData[idx1].push_back (ForceData {force_idx, idx1,
            // idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase,
            // vtype.k});
            // bondTorsionData[idx2].push_back (ForceData {force_idx, idx1,
            // idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase,
            // vtype.k});
            // bondTorsionData[idx3].push_back (ForceData {force_idx, idx1,
            // idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase,
            // vtype.k});
            // bondTorsionData[idx4].push_back (ForceData {force_idx, idx1,
            // idx2, idx3, idx4, 0, 0, vtype.periodicity, vtype.phase,
            // vtype.k});

            //++force_idx;
        }
    }

    // Create the 1-2-3-4 improper terms where 3 is the central atom
    for (auto& dihedral : topology.impropers) {
        dbgmsg("checkpoint10");
        const molib::Atom& atom1 = *get<0>(dihedral);
        const molib::Atom& atom2 = *get<1>(dihedral);
        const molib::Atom& atom3 = *get<2>(dihedral);
        const molib::Atom& atom4 = *get<3>(dihedral);
        const int idx1 = topology.get_index(atom1);
        const int idx2 = topology.get_index(atom2);
        const int idx3 = topology.get_index(atom3);
        const int idx4 = topology.get_index(atom4);
        int type1 = topology.get_type(atom1);
        int type2 = topology.get_type(atom2);
        int type3 = topology.get_type(atom3);
        int type4 = topology.get_type(atom4);

        // Check for modified residues
        int number_protein = (atom1.br().rest() == molib::Residue::protein) +
                             (atom2.br().rest() == molib::Residue::protein) +
                             (atom3.br().rest() == molib::Residue::protein) +
                             (atom4.br().rest() == molib::Residue::protein);

        int number_not_set =
            (atom1.gaff_type() != "???") + (atom2.gaff_type() != "???") +
            (atom3.gaff_type() != "???") + (atom4.gaff_type() != "???");

        // It should all be protein or non-protein, otherwise we fix it
        if ((number_protein == 1 || number_protein == 2 ||
             number_protein == 3) &&
            (number_not_set == 1 || number_not_set == 2 ||
             number_not_set == 3)) {
            retype_amber_protein_atom_to_gaff(atom1, type1);
            retype_amber_protein_atom_to_gaff(atom2, type2);
            retype_amber_protein_atom_to_gaff(atom3, type3);
            retype_amber_protein_atom_to_gaff(atom4, type4);
        }

        try {
            const ForceField::TorsionTypeVec& v_ttype =
                __ffield->get_improper_type(type1, type2, type3,
                                            type4);  // cannot make it const ??

            for (auto& ttype : v_ttype) {
                dbgmsg("force_idx = " << force_idx << " idx1 = " << idx1
                                      << " idx2 = " << idx2 << " idx3 = "
                                      << idx3 << " idx4 = " << idx4);
                bondTorsion->addTorsion(idx1, idx2, idx3, idx4,
                                        ttype.periodicity, ttype.phase,
                                        ttype.k);
                bondTorsionData[idx1].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});
                bondTorsionData[idx2].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});
                bondTorsionData[idx3].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});
                bondTorsionData[idx4].push_back(
                    ForceData{force_idx, idx1, idx2, idx3, idx4, 0, 0,
                              ttype.periodicity, ttype.phase, ttype.k});

                ++force_idx;
            }
        } catch (ParameterError& e) {
            dbgmsg(e.what() << " (WARNINGS ARE NOT INCREASED)");
        }
    }

    if (warn > 0) {
        throw Error("die : missing parameters detected");
    }
}

void SystemTopology::init_positions(const geometry::Point::Vec& crds) {
    dbgmsg("entering init_positions crds.size() = " << crds.size());
    vector<OpenMM::Vec3> positions_in_nm;

    for (auto& crd : crds) {
        positions_in_nm.push_back(OpenMM::Vec3(
            crd.x() * OpenMM::NmPerAngstrom, crd.y() * OpenMM::NmPerAngstrom,
            crd.z() * OpenMM::NmPerAngstrom));
    }

    context->setPositions(positions_in_nm);

    dbgmsg("exiting init_positions");
}

geometry::Point::Vec SystemTopology::get_positions_in_nm() {
    // The true parameter is to enforce Periodic Boundary Conditions
    auto positions_in_nm =
        context->getState(OpenMM::State::Positions).getPositions();
    geometry::Point::Vec result;
    result.reserve(positions_in_nm.size());
    for (size_t i = 0; i < positions_in_nm.size(); ++i) {
        result.emplace_back(geometry::Point(positions_in_nm[i][0],
                                            positions_in_nm[i][1],
                                            positions_in_nm[i][2]));
    }

    return result;
}

geometry::Point::Vec SystemTopology::get_forces() {
    auto forces_in_nm =
        context->getState(OpenMM::State::Forces, true).getForces();
    geometry::Point::Vec result;
    result.reserve(forces_in_nm.size());
    for (size_t i = 0; i < forces_in_nm.size(); ++i) {
        result.emplace_back(geometry::Point(
            forces_in_nm[i][0], forces_in_nm[i][1], forces_in_nm[i][2]));
    }

    return result;
}

double SystemTopology::get_potential_energy() {
    return context->getState(OpenMM::State::Energy, true).getPotentialEnergy();
}

double SystemTopology::get_kinetic_energy() {
    return context->getState(OpenMM::State::Energy, true).getKineticEnergy();
}

void SystemTopology::print_box_vector_size() {
    OpenMM::Vec3 x, y, z;
    system->getDefaultPeriodicBoxVectors(x, y, z);
    std::cerr << "x " << x[0] << " " << x[1] << " " << x[2] << std::endl;
    std::cerr << "y " << y[0] << " " << y[1] << " " << y[2] << std::endl;
    std::cerr << "z " << z[0] << " " << z[1] << " " << z[2] << std::endl;
}

void SystemTopology::minimize(const double tolerance,
                              const int max_iterations) {
    log_note << "Minimizing with " << max_iterations
             << " iterations with a tolerence of " << tolerance << "\n";
    OpenMM::LocalEnergyMinimizer::minimize(*context, tolerance, max_iterations);
}

double SystemTopology::get_energies() {
    return get_potential_energy() + get_kinetic_energy();
}

void SystemTopology::set_temperature() {
    context->setVelocitiesToTemperature(300);
}

void SystemTopology::set_box_vector() {
    context->setPeriodicBoxVectors(OpenMM::Vec3(6, 0, 0), OpenMM::Vec3(0, 6, 0),
                                   OpenMM::Vec3(0, 0, 6));
}

void SystemTopology::print_energies() {
    cerr << "Potential Eneriges: " << get_potential_energy()
         << "\tKinetic Energies: " << get_kinetic_energy()
         << "\tTotal Energies: " << get_energies() << endl;
}

void SystemTopology::dynamics(const int steps) {
    if (__integrator_used == integrator_type::verlet && __thermostat_idx == -1)
        log_warning << "No thermostat set, performing NVE dynamics" << endl;

    auto start = high_resolution_clock::now();

    integrator->step(steps);

    auto length = duration_cast<seconds>(high_resolution_clock::now() - start);
    cerr << "Time taken by function: " << length.count() << " seconds\n";

    print_energies();
    save_checkpoint();
}

void SystemTopology::load_checkpoint(const std :: string & checkpoint) {
    try {
        ifstream file(checkpoint, ios::in | ios::binary);
        if(file.is_open()){
            context->loadCheckpoint(file);
            file.close();
        }
        else {
            cerr << "Error could not open file: " << checkpoint << endl;
            exit(1);
        }
    } catch (const std::exception& e) {
        log_error << "Checkpoint failed to load " << e.what() << endl;
        exit(1);
    }
}

void SystemTopology::save_checkpoint() {
    string checkpoint = "checkpoint" + to_string(checkpoint_num % num_checkpoints) + ".chk";
    cerr << "Writing checkpoint to file " << checkpoint << endl;
    ofstream file(checkpoint, ios::out | ios::binary | ios::trunc);
    context->createCheckpoint(file);
    file.close();
}

void SystemTopology::save_checkpoint_candock(int x, int y)
{
    string checkpoint = "checkpoint" + to_string(checkpoint_num % num_checkpoints) + ".chk_candock";
    ofstream file(checkpoint, ios::out | ios::trunc);
    file << x;
    file << "\n";
    file << y;
    checkpoint_num++;
    file.close();
}


}  // namespace OMMIface
}  // namespace statchem
