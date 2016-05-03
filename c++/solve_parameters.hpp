#pragma once

namespace realevol {

enum ode_solve_method {RungeKutta, Lanczos};
using operator_t = triqs::operators::many_body_operator_generic<time_expr>;

// All the arguments of the solve function
struct solve_parameters_t {

 /// Hamiltonian
 operator_t h;

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 int verbosity = ((triqs::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Planck constant
 double hbar = 1.0;

 /// Method to solve the Schroedinger equation
 ode_solve_method method = Lanczos;

 solve_parameters_t() {}
 solve_parameters_t(operator_t const& h) : h(h) {}
};

}
