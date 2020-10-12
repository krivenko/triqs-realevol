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

#include <triqs/utility/variant.hpp>

#include "common.hpp"
#include "init_state.hpp"
#include "propagator.hpp"

namespace realevol {

// All the arguments of the compute_2t_obs function
struct compute_2t_obs_parameters_t {

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 int verbosity = ((triqs::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Compute lesser Green's function
 bool compute_g_l = true;

 /// Compute greater Green's function
 bool compute_g_g = true;

 /// Compute susceptibility
 bool compute_chi = true;

 /// Compute components of observables with the first time argument within this range
 std::pair<double, double> t_range = {-INFINITY, INFINITY};

 /// Compute components of observables with the second time argument within this range
 std::pair<double, double> tp_range = {-INFINITY, INFINITY};

 /// Compute components of observables with time arguments satisfying |t-t'|<=delta_t_max
 double delta_t_max = INFINITY;

 /// Planck's constant
 double hbar = 1.0;

 /// Hamiltonian interpolation between time slices
 h_interpolation hamiltonian_interpol = Rectangle;

 /// Use Lanczos algorithm to exponentiate matrices of this size or bigger
 int lanczos_min_matrix_size = 11;

 /// Lanczos convergence threshold for the GS energy, for each invariant subspace
 std::map<long,double> lanczos_gs_energy_tol = std::map<long,double>({});

 /// Maximal dimension of the Krylov space, for each invariant subspace
 std::map<long,int> lanczos_max_krylov_dim = std::map<long,int>({});

 compute_2t_obs_parameters_t() = default;
};

}
