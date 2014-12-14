#pragma once

namespace realevol {

enum ode_solve_method {method_runge_kutta, method_lanczos};

// All the arguments of the solve function
template<bool ComplexOp> struct solve_parameters_t {

 /// Hamiltonian
 operator_t<ComplexOp> h;

 /// Verbosity level
 int verbosity = ((boost::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Observables to be measured
 std::map<std::string,operator_t<ComplexOp>> observables = {};

 /// Planck constant
 double hbar = 1.0;

 /// Time mesh to solve the Schroedinger equation
 any_mesh_t mesh;

 /// Method to solve the Schroedinger equation
 ode_solve_method method = method_lanczos;

 /// Keep in memory the state vector on this number of time points
 long stored_psi_values = 10;

 /// Mesh downsampling factors for the observables
 //std::map<std::string,int> mesh_downsampling = (std::map<std::string,int>{});

 solve_parameters_t(operator_t<ComplexOp> const& h, any_mesh_t const& mesh) : h(h), mesh(mesh) {}
};
}
