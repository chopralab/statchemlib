#ifndef _STCH_ALLSCORE_POSE_HPP_
#define _STCH_ALLSCORE_POSE_HPP_

#include "Program.hpp"

namespace statchem_prog {

class AllScorePose : public Program {
   public:
    AllScorePose();

    virtual int run(int argc, char* argv[]);
};


template<> ProgramInfo program_information<AllScorePose>();

}

#endif
