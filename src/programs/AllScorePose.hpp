#ifndef _STCH_ALLSCORE_POSE_HPP_
#define _STCH_ALLSCORE_POSE_HPP_

#include "Program.hpp"

#include <string>
#include "statchem/molib/molecules.hpp"

namespace statchem_prog {

class AllScorePose : public Program {
   public:
    AllScorePose();

    virtual bool process_options(int argc, char* argv[]) override;
    virtual int run() override;
   private:
    std::string __dist;
    statchem::molib::Molecules __receptor_mols;
    statchem::molib::Molecules __ligand_mols;
    size_t __num_threads;
    bool __constant_receptor;
};


template<> ProgramInfo program_information<AllScorePose>();

}

#endif
