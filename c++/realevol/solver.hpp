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

#include <string>
#include <vector>
#include <map>
#include <utility>

#include <triqs/mpi/base.hpp>

#include "triqs/operators/many_body_operator.hpp"
#include "compute_2t_obs_parameters.hpp"

namespace realevol {

using namespace triqs::gfs;

using indices_type = operators::indices_t;
using chi_indices_t = std::vector<std::pair<std::string,utility::variant_int_string>>;

class solver {

 gf_struct_t gf_struct;                              // Block structure of the Green functions
 chi_indices_t chi_indices;                          // Indices of susceptibilities
 block_gf_2t_t g_l, g_g;                             // Lesser and greater GF's to be calculated
 gf_2t_t chi;                                        // Susceptibility to be calculated
 init_state const * initial_state = nullptr;         // Initial state at t=0
 triqs::mpi::communicator comm;                      // MPI communicator
 compute_2t_obs_parameters_t compute_2t_obs_params;  // Parameters of the last call to solve
 gf_mesh<retime> t_mesh;                             // 1D time mesh to use in calculations

public:

 solver(gf_struct_t const& gf_struct, chi_indices_t const& chi_indices, double t_max, int n_t);

 /// Compute observables that are functions of two times
 template<typename HamiltonianType>
 void compute_2t_obs(HamiltonianType const& h, compute_2t_obs_parameters_t const& p);

 /// Get the initial state at t=0
 init_state const& get_initial_state() const {
  if(initial_state == nullptr) TRIQS_RUNTIME_ERROR << "Initial state is not set!";
  return *initial_state;
 }

 /// Set the initial state
 void set_initial_state(init_state const& initial_state) { this->initial_state = &initial_state; }

 /// Set of parameters used in the last call to compute_2t_obs()
 compute_2t_obs_parameters_t get_last_compute_2t_obs_parameters() const { return compute_2t_obs_params; }

 /// Lesser GF in real time
 block_gf_2t_view get_g_l() { return g_l; }

 /// Greater GF in real time
 block_gf_2t_view get_g_g() { return g_g; }

 /// Susceptibility in real time
 gf_2t_view get_chi() { return chi; }

};

}
