#include <statchem/capi/capi.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char** argv) {
    clock_t start = clock();

    if (!initialize_complex("files/1aaq.pdb")) {
        printf("%s\n", cd_get_error());
        return 1;
    }

    size_t rec_atom_count = receptor_atom_count();
    if (rec_atom_count != 1516) {
        return 2;
    }

    char* chain = (char*)malloc(rec_atom_count * sizeof(char));
    size_t* resi = (size_t*)malloc(rec_atom_count * sizeof(size_t));
    size_t* rest = (size_t*)malloc(rec_atom_count * sizeof(size_t));
    char* resn = (char*)malloc(rec_atom_count * sizeof(char));
    size_t* elem = (size_t*)malloc(rec_atom_count * sizeof(size_t));

    if (!receptor_atom_details(chain, resi, rest, resn, elem)) {
        printf("%s\n", cd_get_error());
        return 3;
    }

    size_t lig_atom_count = ligand_atom_count();
    if (lig_atom_count != 42) {
        return 4;
    }

    char* lchain = (char*)malloc(lig_atom_count * sizeof(char));
    size_t* lresi = (size_t*)malloc(lig_atom_count * sizeof(size_t));
    size_t* lrest = (size_t*)malloc(lig_atom_count * sizeof(size_t));
    size_t* lelem = (size_t*)malloc(lig_atom_count * sizeof(size_t));
    int* ltype = (size_t*)malloc(lig_atom_count * sizeof(int));
    if (!ligand_atom_details(lchain, lresi, lrest, lelem, ltype)) {
        printf("%s\n", cd_get_error());
        return 5;
    }

    free(chain);
    free(resi);
    free(resn);
    free(rest);
    free(elem);
    free(lchain);
    free(lresi);
    free(lrest);
    free(lelem);
    free(ltype);

    return 0;
}
