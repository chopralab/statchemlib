#include <statchem/capi/capi.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char** argv) {

    if (!initialize_complex("files/benzene.mol2")) {
        printf("%s\n", cd_get_error());
        return 1;
    }

    ligand_ligand_atoms();

    size_t lig_atom_count = ligand_atom_count();
    if (lig_atom_count != 6) {
        return 2;
    }

    char* lchain = (char*)malloc(lig_atom_count * sizeof(char));
    size_t* lresi = (size_t*)malloc(lig_atom_count * sizeof(size_t));
    size_t* lrest = (size_t*)malloc(lig_atom_count * sizeof(size_t));
    size_t* lelem = (size_t*)malloc(lig_atom_count * sizeof(size_t));
    int* ltype = (size_t*)malloc(lig_atom_count * sizeof(int));

    if (!ligand_atom_details(lchain, lresi, lrest, lelem, ltype)) {
        printf("%s\n", cd_get_error());
        return 3;
    }

    for (size_t i = 0; i < lig_atom_count; ++i) {
        if (ltype[i] != 22) {
            return 4;
        }
    }

    size_t neighbors[2];
    ligand_get_neighbors(0, neighbors);
    neighbors[0] = neighbors[1] = 9;
    ligand_get_neighbors(1, neighbors);

    // Change the connectivity
    if (!remove_ligand_bond(1, 2)) {
        return 5;
    }

    if (is_adjacent(1, 2)) {
        return 10;
    }

    if (ligand_bond_count() != 5) {
        return 7;
    }

    neighbors[0] = neighbors[1] = 9;
    ligand_get_neighbors(0, neighbors);

    neighbors[0] = neighbors[1] = 9;
    ligand_get_neighbors(1, neighbors);

    ligand_ligand_atoms();

    if (!ligand_atom_details(lchain, lresi, lrest, lelem, ltype)) {
        printf("%s\n", cd_get_error());
        return 3;
    }

    for (size_t i = 0; i < lig_atom_count; ++i) {
        if (ltype[i] == 22) {
            return 6;
        }
    }

    return 0;
}
