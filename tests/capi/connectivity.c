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

    if (!is_adjacent(1519, 1520)) {
        return 1;
    }

    if (is_adjacent(1519, 1521)) {
        return 1;
    }

    if (ligand_bond_count() != 41) {
        return 2;
    }

    // Change the connectivity
    if (!remove_ligand_bond(1519, 1520)) {
        return 1;
    }

    if (ligand_bond_count() != 40) {
        return 3;
    }

    if (!add_ligand_bond(1519, 1521)) {
        return 1;
    }

    // Logic is flipped
    if (is_adjacent(1519, 1520)) {
        return 1;
    }

    if (!is_adjacent(1519, 1521)) {
        return 1;
    }

    return 0;
}
