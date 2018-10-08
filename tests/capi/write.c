#include <statchem/capi/capi.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char** argv) {

    if (!initialize_complex("files/1aaq.pdb")) {
        printf("%s\n", cd_get_error());
        return 1;
    }

    if (!write_complex("temp.pdb")) {
        printf("%s\n", cd_get_error());
        return 1;
    }

    if (!initialize_complex("temp.pdb")) {
        printf("%s\n", cd_get_error());
        return 1;
    }

    remove("temp.pdb");

    size_t rec_atom_count = receptor_atom_count();
    if (rec_atom_count != 1516) {
        return 2;
    }

    size_t lig_atom_count = ligand_atom_count();
    if (lig_atom_count != 42) {
        return 2;
    }

    size_t lig_bond_count = ligand_bond_count();
    if (lig_bond_count != 41) {
        return 3;
    }

    return 0;
}
