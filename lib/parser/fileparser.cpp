/* This is fileparser.cpp and is part of StatChemLIB
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

#include "statchem/parser/fileparser.hpp"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "statchem/fragmenter/fragmenter.hpp"
#include "statchem/geometry/coordinate.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/molib/bond.hpp"
#include "statchem/molib/nrset.hpp"

using namespace std;
using namespace statchem::molib;

namespace statchem {
namespace parser {
void FileParser::set_flags(unsigned int hm) { p->set_hm(hm); }
bool FileParser::parse_molecule(Molecules& mols) {
    p->parse_molecule(mols);
    dbgmsg("PARSED MOLECULES : " << endl << mols);
    return !mols.empty();
}

Molecules FileParser::parse_molecule() {
    Molecules mols;
    p->parse_molecule(mols);
    dbgmsg("PARSED MOLECULES : " << endl << mols);

    if (mols.empty()) {
        throw Error("die : could not read molecules from file");
    }

    return mols;
}

FileParser::FileParser(const string& molecule_file, unsigned int hm,
                       const int num_occur) {
    prepare_parser(molecule_file, hm, num_occur);
}

void FileParser::prepare_parser(const string& molecule_file, unsigned int hm,
                                const int num_occur) {
    if (fileio::file_size(molecule_file) <= 0) {
        throw Error("die : file not valid: " + molecule_file +
                    ". Check to see if it exists and has contents!");
    }

    auto ret = molecule_file.find_last_of(".");

    if (ret == string::npos) {
        throw Error(
            "die : could not determine the file type of the input molecule");
    }

    string extension = molecule_file.substr(ret + 1);
    transform(extension.begin(), extension.end(), extension.begin(), ::toupper);

    std::shared_ptr<istream> temp_molecule_stream =
        std::make_shared<ifstream>(molecule_file, std::ios::in);

    prepare_parser(temp_molecule_stream, extension, hm, num_occur);
}

void FileParser::prepare_parser(std::shared_ptr<std::istream>& stream,
                                const std::string& extension, unsigned int hm,
                                const int num_occur) {
    molecule_stream = stream;

    dbgmsg("pdb reader options = "
           << hm << " molecule is a " << boolalpha
           << (extension == "PDB" || extension == "ENT"
                   ? "PDB file"
                   : (extension == "MOL2" ? "Mol2 file"
                                          : "undetermined file type")));

    if (extension == "PDB" || extension == "ENT") {
        p = std::unique_ptr<Parser>(
            new PdbParser(*molecule_stream, hm, num_occur));
    } else if (extension == "MOL2") {
        p = std::unique_ptr<Parser>(
            new Mol2Parser(*molecule_stream, hm, num_occur));
    } else {
        throw Error(
            "die : could not determine the file type of the input molecule");
    }
}
}
}
