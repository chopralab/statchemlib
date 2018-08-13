#ifndef _STCH_PHYS_MINIMIZE_HPP_
#define _STCH_PHYS_MINIMIZE_HPP_

#include "Program.hpp"

#include "statchem/molib/molecules.hpp"
#include "statchem/modeler/forcefield.hpp"

namespace statchem_prog {

class PhysMinimize : public Program {
   public:
    PhysMinimize();

    virtual bool process_options(int argc, char* argv[]) override;
    virtual int run() override;
   private:
    statchem::molib::Molecules __receptor_mols;
    statchem::molib::Molecules __ligand_mols;
    bool __constant_receptor;

    statchem::OMMIface::ForceField __ffield;
    double __mini_tol;
    int __iter_max;
    std::string __platform, __precision, __accelerators;
};


template<> ProgramInfo program_information<PhysMinimize>();

}

#endif
