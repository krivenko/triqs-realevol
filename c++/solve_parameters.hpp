#pragma once

namespace realevol {

// All the arguments of the solve function
template<typename OperatorType, typename Mesh>
struct solve_parameters_t {

 /// Hamiltonian
 OperatorType h;

 /// Verbosity level
 int verbosity = ((boost::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Observables to be measured
 std::map<std::string,OperatorType> observables = {};

 /// Planck constant
 double planck_constant = 1.0;

 // Time mesh to solve the Schroedinger equation
 Mesh mesh;

 //.add_field("ode_solve_method", triqs::params::no_default<ode_solve_method>(), "Method to solve Schroedinger equation");
 //.add_field("mesh_downsampling", dict_t<int>({}), "Mesh downsampling factors for the observables");

 solve_parameters_t(OperatorType const& h, Mesh const& mesh) : h(h), mesh(mesh) {}
};
}
