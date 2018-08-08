#ifndef _STCH_PROGRAM_MANAGER_HPP_
#define _STCH_PROGRAM_MANAGER_HPP_

#include <functional>
#include <memory>
#include <vector>
#include "Program.hpp"

namespace statchem_prog {

using program_creator_t = std::function<std::unique_ptr<Program>()>;

class ProgramManager final {
   private:
    ProgramManager();

   public:
    static ProgramManager& get();

    template <class Program>
    void add_format() {
        auto info = program_information<Program>();

        register_program(info, []() {
            return std::unique_ptr<Program>(
                new Program());  // NOLINT no make_unique in C++11
        });
    }

    program_creator_t name(const std::string& name);

    std::vector<ProgramInfo> programs();

   private:
    using program_map_t =
        std::vector<std::pair<ProgramInfo, program_creator_t>>;
    using iterator = program_map_t::const_iterator;

    static iterator find_name(const program_map_t& programs,
                              const std::string& name);

    void register_program(ProgramInfo info, program_creator_t creator);

    program_map_t __programs;
};
}

#endif
