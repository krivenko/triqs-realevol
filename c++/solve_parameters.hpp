#pragma once

namespace realevol {

enum ode_solve_method {RungeKutta, Lanczos};
using operator_t = realevol::operators::many_body_operator_generic<time_expr>;

// All the arguments of the solve function
struct solve_parameters_t {

 /// Hamiltonian
 operator_t h;

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 int verbosity = ((triqs::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Planck constant
 double hbar = 1.0;

 /// Use a thermal equilibrium state with inverse temperature :math:`\\beta` as the initial state.
 bool thermal_init_state = true;

 /// Operator to generate the initial state
 // :math:`\\hat\\rho_0\\propto\\exp(-\\beta\\hat h_0)`, if `thermal_init_state = True`,
 // :math:`\\psi_0\\rangle = \\hat h_0|vac\\rangle` otherwise
 operator_t h0;

 /// Inverse temperature
 /// default: +inf
 double beta = HUGE_VAL;

 /// Number of binary digits per bosonic degree of freedom
 /// type: dict(Operator index : int)
 std::map<realevol::operators::indices_t, int> bits_per_boson = {};

 /// Method to solve the Schroedinger equation
 ode_solve_method method = Lanczos;

 solve_parameters_t() {}
 solve_parameters_t(operator_t const& h) : h(h) {}
};

}
