#ifndef _STCH_MAKE_OBJECTIVE_HPP_
#define _STCH_MAKE_OBJECTIVE_HPP_

#include "Program.hpp"

#include <memory>
#include <set>

#include "statchem/molib/molecules.hpp"
#include "statchem/score/kbff.hpp"

namespace statchem_prog {

class MakeObjective : public Program {
   public:
    MakeObjective();

    virtual bool process_options(int argc, char* argv[]) override;
    virtual int run() override;
   private:
    std::string __dist;
    std::string __ref;
    std::string __comp;
    std::string __func;
    double __cutoff;
    double __step_size;
    std::unique_ptr<statchem::score::KBFF> __score;

    std::set<int> __atom_types;
    std::string __obj_dir;
};


template<> ProgramInfo program_information<MakeObjective>();

}

#endif
