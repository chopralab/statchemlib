#ifndef _STCH_PHYS_DYNAMICS_HPP_
#define _STCH_PHYS_DYNAMICS_HPP_

#include "Program.hpp"

#include <memory>

#include "statchem/molib/molecules.hpp"
#include "statchem/modeler/forcefield.hpp"

namespace statchem_prog {

class PhysDynamics : public Program {
   public:
    PhysDynamics();

    virtual bool process_options(int argc, char* argv[]) override;
    virtual int run() override;
   private:
    std::string __dist;
    statchem::molib::Molecules __receptor_mols;
    statchem::molib::Molecules __ligand_mols;
    bool __constant_receptor;
    size_t __num_threads;

    statchem::OMMIface::ForceField __ffield;
    std::string __integrator;
    double __mini_tol;
    int __iter_max;
    double __dist_cut;
    double __temperature;
    double __friction;
    int __dynamics_steps;
    double __dynamics_step_size;
    std::string __platform, __precision, __accelerators, __checkpoint;
};


template<> ProgramInfo program_information<PhysDynamics>();

}

#endif
