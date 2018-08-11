#ifndef _STCH_KB_MINIMIZE_HPP_
#define _STCH_KB_MINIMIZE_HPP_

#include "Program.hpp"

namespace statchem_prog {

class KBMinimize : public Program {
   public:
    KBMinimize();

    virtual int run(int argc, char* argv[]);
};


template<> ProgramInfo program_information<KBMinimize>();

}

#endif
