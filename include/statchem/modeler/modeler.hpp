#ifndef MODELER_H
#define MODELER_H
#include <map>
#include <set>
#include <string>
#include <vector>

#include "statchem/geometry/coordinate.hpp"
#include "statchem/modeler/systemtopology.hpp"
#include "statchem/modeler/topology.hpp"
#include "statchem/molib/molecule.hpp"

namespace statchem {

namespace score {
class Score;
}

namespace OMMIface {

class Modeler {
   public:
    class MinimizationError : public Error {
       public:
        MinimizationError(const std::string& msg) : Error(msg) {}
    };

   private:
    const ForceField* __ffield;
    std::string __fftype;
    double __scale;
    double __tolerance;
    int __max_iterations;
    bool __use_constraints;
    double __step_size_in_ps;
    double __temperature;
    double __friction;

    geometry::Point::Vec __positions;
    Topology __topology;

    int __dynamics_steps;

   public:
    SystemTopology __system_topology;

    Modeler(const ForceField& ffield, const std::string& fftype = "none",
            double scale = 1.0, double tolerance = 0.0001,
            int max_iterations = 100, bool use_constraints = false,
            double step_size_in_fs = 2.0, double temperature = 300.0,
            double friction = 91.0);

    void mask(const molib::Atom::Vec& atoms);
    void unmask(const molib::Atom::Vec& atoms);

    void add_topology(const molib::Atom::Vec& atoms);
    void add_crds(const molib::Atom::Vec& atoms,
                  const geometry::Point::Vec& crds);
    void add_random_crds(const molib::Atom::Vec& atoms);

    geometry::Point::Vec get_state(const molib::Atom::Vec& atoms);

#ifndef NDEBUG
    void minimize_knowledge_based(molib::Molecule& ligand,
                                  molib::Molecule& receptor,
                                  score::Score& score);
#endif

    void minimize_state();
    void dynamics();

    void init_openmm_positions();
    void init_openmm(const std::string& platform,
                     const std::string& precision = "double",
                     const std::string& accelerators = "",
                     SystemTopology::integrator_type type =
                         SystemTopology::integrator_type::none);

    void set_max_iterations(const int max_iterations) {
        __max_iterations = max_iterations;
    }
    void set_num_steps_to_run(const int num_steps_to_run) {
        __dynamics_steps = num_steps_to_run;
    }

    double potential_energy();

    const ForceField& get_forcefield() { return *__ffield; }
};
}  // namespace OMMIface
}  // namespace statchem

#endif
