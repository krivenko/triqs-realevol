/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include "common.hpp"

namespace realevol {

enum ode_solve_method {RungeKutta, Lanczos};

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
