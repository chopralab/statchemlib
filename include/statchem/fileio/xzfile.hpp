#ifndef XZFILE_H
#define XZFILE_H

#include <sstream>

namespace statchem {
namespace fileio {

std::stringstream read_xzfile(const std::string& filename);

}
}

#endif
