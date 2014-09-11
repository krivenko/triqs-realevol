#include "realevol.hpp"
#include "uniform_mesh.hpp"

namespace realevol {

#define SOLVER_METHOD(M)    template<typename Mesh, bool ComplexOperators> auto solver<Mesh,ComplexOperators>::M 

SOLVER_METHOD(solve(operator_t h, parameters_t params, dict_t<operator_t> observables) -> void){
}

SOLVER_METHOD(solve_parameters() -> parameters_t) {
    boost::mpi::communicator world;
    auto pdef = parameters_t{};

    pdef.add_field("verbosity", (world.rank()==0 ? int(3) : int(0)), "Verbosity level")
        //.add_field("mesh", triqs::params::no_default<Mesh>(), "Time mesh to solve the Schroedinger equation")
        .add_field("planck_constant", double(1.0), "Planck constant");
        //.add_field("mesh_downsampling", dict_t<int>({}), "Mesh downsampling factors for the observables");

    return pdef;
}

SOLVER_METHOD(help() -> void) {
// TODO
}

// Explicit instantiations
// Add more mesh types later if needed
template class solver<uniform_mesh<>,false>;
template class solver<uniform_mesh<>,true>;

}