#ifndef _STCH_PROGRAM_HPP_
#define _STCH_PROGRAM_HPP_

#include <stdexcept>
#include <string>
#include <sstream>

namespace statchem_prog {

class Program {
   public:
    virtual bool process_options(int argc, char* argv[]) = 0;
    virtual int run() = 0;

    std::string get_help() const {
        return __help_text.str();
    }

   protected:
    std::stringstream __help_text;
};

class ProgramInfo {
   public:
    ProgramInfo(std::string name) : __name(std::move(name)) {
        if (__name == "") {
            throw std::length_error("a format name can not be an empty string");
        }
    }

    const std::string& name() const { return __name; }

    ProgramInfo& description(std::string description) {
        __description = std::move(description);
        return *this;
    }

    const std::string& description() const { return __description; }

   private:
    std::string __name;
    std::string __description;
};

template <class Program>
ProgramInfo program_information() {
    throw std::logic_error(
        "program_information is unimplemented for this format");
}
}

#endif
