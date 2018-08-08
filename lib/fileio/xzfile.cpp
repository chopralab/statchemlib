#include "statchem/fileio/xzfile.hpp"
#include "statchem/helper/error.hpp"

#include <lzma.h>
#include <fstream>
#include <vector>

namespace statchem {
namespace fileio {

std::stringstream read_xzfile(const std::string& filename) {
    
    std::ifstream checking(filename, std::ifstream::binary);
    std::vector<char> src((std::istreambuf_iterator<char>(checking)),
                           std::istreambuf_iterator<char>());
    
    std::vector<char> dst(src.size());

    size_t dataLength = src.size();
    size_t expdLength = dataLength * 5;

    bool done = false;

    lzma_stream strm = LZMA_STREAM_INIT;
    strm.next_in = reinterpret_cast<unsigned char*>(src.data());
    strm.avail_in = dataLength;
    strm.total_out = 0;

    lzma_ret status = lzma_stream_decoder(&strm, UINT64_MAX, 0);

    if (status != LZMA_OK) {
	    throw Error("Problem with Initializing LZMA");
    }

    lzma_action action = LZMA_RUN;

    do {

	    // extend decompressed if too short
	    if (strm.total_out >= dst.size()) {
		    dst.resize(dst.size() + expdLength);
	    }

	    strm.next_out = reinterpret_cast<unsigned char*>(dst.data() + strm.total_out);
    	strm.avail_out = dst.size() - strm.total_out;

        status = lzma_code(&strm, action);

        if (status == LZMA_STREAM_END) {
		    done = true;
        } else if (status != LZMA_OK) {
		    lzma_end(&strm);
            throw Error("Problem with decompressing LZMA");
	    }
    } while(!done);

    lzma_end (&strm);

    // set actual length, shrinking the vector is
    // very efficient in C++11!
    dst.resize(strm.total_out);

    std::stringstream out(std::string(dst.data(), dst.size()));

    return out;
}

}
}
