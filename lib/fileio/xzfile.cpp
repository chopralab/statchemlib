/* This is xzfile.cpp and is part of StatChemLIB
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "statchem/fileio/xzfile.hpp"
#include "statchem/helper/error.hpp"

#include <lzma.h>
#include <fstream>
#include <vector>

namespace statchem {
namespace fileio {

std::string read_xzfile(const std::string& filename) {
    
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

    return std::string(dst.data(), dst.size());
}

}
}
