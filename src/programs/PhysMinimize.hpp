#ifndef _STCH_PHYS_MINIMIZE_HPP_
#define _STCH_PHYS_MINIMIZE_HPP_

#include "Program.hpp"

namespace statchem_prog {

class PhysMinimize : public Program {
   public:
    PhysMinimize();

    virtual int run(int argc, char* argv[]);
};


template<> ProgramInfo program_information<PhysMinimize>();

}

#endif
