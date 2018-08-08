#ifndef _STCH_SCORE_POSE_HPP_
#define _STCH_SCORE_POSE_HPP_

#include "Program.hpp"

namespace statchem_prog {

class ScorePose : public Program {
   public:
    ScorePose();

    virtual int run(int argc, char* argv[]);
};


template<> ProgramInfo program_information<ScorePose>();

}

#endif
