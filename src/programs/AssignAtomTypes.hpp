#ifndef _STCH_ASSIGNATOMTYPES_POSE_HPP_
#define _STCH_ASSIGNATOMTYPES_POSE_HPP_

#include "Program.hpp"

#include <string>
#include "statchem/molib/molecules.hpp"

namespace statchem_prog {

class AssignAtomTypes : public Program {
   public:
    AssignAtomTypes();

    virtual bool process_options(int argc, char* argv[]) override;
    virtual int run() override;
   private:
    statchem::molib::Molecules __molecules;
    std::string __output;
};


template<> ProgramInfo program_information<AssignAtomTypes>();

}

#endif
