/* This is parser.hpp and is part of StatChemLIB
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

#ifndef PARSER_H
#define PARSER_H

#include <string>
#include "statchem/molib/molecules.hpp"

namespace statchem {

namespace parser {
enum pdb_read_options {
    first_model = 1 << 0,
    all_models = 1 << 1,
    hydrogens = 1 << 2,
    skip_hetatm = 1 << 3,
    skip_atom = 1 << 4,
    sparse_macromol = 1 << 5,
    docked_poses_only = 1 << 6,
    protein_poses_only = 1 << 7,
};

class Parser {
   protected:
    std::istream& __stream;
    std::mutex __concurrent_read_mtx;
    unsigned int __hm;
    const int __num_occur;
    bool __giant_molecule;
    void __generate_molecule(molib::Molecules&, bool&, const std::string&);
    void __generate_assembly(molib::Molecules&, bool&, int, const std::string&);
    void __generate_model(molib::Molecules&, bool&, int);

   public:
    Parser(std::istream& molecule_file, unsigned int hm = all_models,
           const int num_occur = -1)
        : __stream(molecule_file),
          __hm(hm),
          __num_occur(num_occur),
          __giant_molecule(false) {}
    virtual ~Parser() {}
    virtual void parse_molecule(molib::Molecules&) = 0;
    virtual void set_pos(std::streampos pos);
    virtual void set_hm(unsigned int hm);
};
}
}

#endif  // PARSER_H
