#ifndef PATH_H
#define PATH_H

#include <string>

namespace statchem {

class Path {
   public:
    static std::string join(const std::string& str1, const std::string& str2);
};
}

#endif
