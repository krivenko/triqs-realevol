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
#include "init_state.hpp"

namespace realevol {

// All the arguments of the compute_2t_obs function
struct compute_2t_obs_parameters_t {

 /// Hamiltonian
 operator_t h;

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 int verbosity = ((triqs::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Planck constant
 double hbar = 1.0;

 /// Compute retarded and advanced Green's functions
 bool compute_gf_ret_adv = true;

 compute_2t_obs_parameters_t() = default;
 compute_2t_obs_parameters_t(operator_t const& h) : h(h) {}
};

}
