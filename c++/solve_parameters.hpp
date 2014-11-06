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

 /// Time mesh to solve the Schroedinger equation
 Mesh mesh;

 /// Method to solve the Schroedinger equation
 enum ode_solve_method {runge_kutta, lanzcos} method = lanzcos;

 /// Mesh downsampling factors for the observables
 std::map<std::string,int> mesh_downsampling = (std::map<std::string,int>{});

 solve_parameters_t(OperatorType const& h, Mesh const& mesh) : h(h), mesh(mesh) {}
};
}
