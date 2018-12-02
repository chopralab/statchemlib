#include "statchem/fileio/fileout.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/molib/bond.hpp"

namespace statchem {
namespace fileio {

void print_complex_pdb(std::ostream& ss, const molib::Molecule& ligand,
                       const molib::Molecule& receptor, const double energy,
                       const double potential, const int model,
                       const size_t max_clq_id, const double rmsd) {
    ss << "MODEL    " << model << "\n";

    ss << "REMARK   1 MINIMIZED COMPLEX OF " << ligand.name() << " AND "
       << receptor.name() << " WITH SCORE OF " << std::setprecision(6) << energy
       << std::endl;

    if (!std::isnan(potential))
        ss << "REMARK   2 POTENTIAL ENERGY OF " << ligand.name() << " IS "
           << potential << "\n";

    if (max_clq_id != 0xFFFFFF)
        ss << "REMARK   3 THIS IS CONFIGUATION " << ligand.name() << " NUMBER "
           << max_clq_id << "\n";

    if (!std::isnan(rmsd))
        ss << "REMARK   3 DOCKED CONFORMATION HAS AN RMSD OF " << rmsd
           << " FROM ORIGINAL STRUCTURE"
           << "\n";

    int reenum = 0;

    for (auto& patom : receptor.get_atoms()) {
        ss << std::setw(66) << std::left << *patom;
        reenum = patom->atom_number();
    }

    for (auto& presd : receptor.get_residues()) {
        if (help::standard_residues.count(presd->resn()) != 0) continue;

        // Print internal bonds
        for (auto pbond : get_bonds_in(presd->get_atoms())) {
            ss << pbond->get_pdb_remark();
        }

        // Print external bonds
        for (auto pbond : get_bonds_in(presd->get_atoms(), false)) {
            ss << pbond->get_pdb_remark();
        }
    }

    ss << "TER"
       << "\n";

    for (auto& patom : ligand.get_atoms()) {
        patom->set_atom_number(++reenum);
        ss << std::setw(66) << std::left << *patom;
    }

    for (auto& pbond : get_bonds_in(ligand.get_atoms())) {
        ss << pbond->get_pdb_remark();
    }

    ss << "ENDMDL"
       << "\n";
}
}  // namespace fileio
}  // namespace statchem
