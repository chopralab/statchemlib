#ifndef _STCH_SCORE_POSE_HPP_
#define _STCH_SCORE_POSE_HPP_

#include "Program.hpp"

#include <memory>

#include "statchem/molib/molecules.hpp"
#include "statchem/score/score.hpp"

namespace statchem_prog {

class ScorePose : public Program {
   public:
    ScorePose();

    virtual bool process_options(int argc, char* argv[]) override;
    virtual int run() override;
   private:
    std::string __dist;
    statchem::molib::Molecules __receptor_mols;
    statchem::molib::Molecules __ligand_mols;
    bool __constant_receptor;
    size_t __num_threads;

    std::string __ref;
    std::string __comp;
    std::string __func;
    double __cutoff;
    std::unique_ptr<statchem::score::Score> __score;

};


template<> ProgramInfo program_information<ScorePose>();

}

#endif
