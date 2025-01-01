/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2025, I. Krivenko, M. Danilov, P. Kubiczek
 *
 * realevol is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * realevol is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * realevol. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <array>
#include <cstdlib>
#include <cmath>
#include <map>
#include <utility>

#include <mpi/mpi.hpp>

#include "array_utility.hpp"
#include "propagator.hpp"

namespace realevol {

// Parameters of the short time propagation solver
template<std::size_t NPoints> struct solver_parameters_t {

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 int verbosity = ((mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Planck's constant
 double hbar = 1.0;

 /// Hamiltonian interpolation between time slices
 h_interpolation hamiltonian_interpol = Rectangle;

 /// Time argument range restrictions, one per correlator point
 std::array<std::pair<double, double>, NPoints> t_ranges =
    make_array_repeat<std::pair<double, double>, NPoints>(std::make_pair(-INFINITY, INFINITY));

 /// Maximal separations of successive time arguments
 std::array<double, NPoints - 1> delta_t_max = make_array_repeat<double, NPoints - 1>(INFINITY);

 /// Use Lanczos algorithm to exponentiate matrices of this size or bigger
 int lanczos_min_matrix_size = 11;

 /// Lanczos convergence threshold for the GS energy, for each invariant subspace
 std::map<long, double> lanczos_gs_energy_tol = {};

 /// Maximal dimension of the Krylov space, for each invariant subspace
 std::map<long, int> lanczos_max_krylov_dim = {};

 solver_parameters_t() = default;
};

}
