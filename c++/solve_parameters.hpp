#pragma once

namespace realevol {

template<bool ComplexOperators>
using operator_coeff_t = typename std::conditional<ComplexOperators,time_expr_c,time_expr_r>::type;
template<bool ComplexOperators>
using operator_t = triqs::utility::many_body_operator<operator_coeff_t<ComplexOperators>>;

enum ode_solve_method {runge_kutta_method, lanczos_method};

// All the arguments of the solve function
template<bool ComplexOperators>
struct solve_parameters_t {

 /// Hamiltonian
 operator_t<ComplexOperators> h;

 /// Verbosity level
 int verbosity = ((boost::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Observables to be measured
 std::map<std::string,operator_t<ComplexOperators>> observables = {};

 /// Planck constant
 double planck_constant = 1.0;

 /// Time mesh to solve the Schroedinger equation
 any_mesh_t mesh;

 /// Method to solve the Schroedinger equation
 ode_solve_method method = lanczos_method;

 /// Keep in memory the state vector on this number of time points
 long psi_time_points = 10;

 /// Mesh downsampling factors for the observables
 std::map<std::string,int> mesh_downsampling = (std::map<std::string,int>{});

 solve_parameters_t(operator_t<ComplexOperators> const& h, any_mesh_t const& mesh) : h(h), mesh(mesh) {}
};
}
