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
#include <triqs/gfs.hpp>

#include "time_expr.hpp"
#include "triqs/operators/many_body_operator.hpp"
#include "compute_gf_parameters.hpp"

namespace realevol {

using namespace triqs::gfs;

using indices_type = operators::indices_t;
using gf_2t_t = gf<cartesian_product<retime, retime>>;
using block_gf_2t_t = block_gf<cartesian_product<retime,retime>>;

class solver {

 gf_struct_t gf_struct, chi_struct;          // Block structure of the Green functions and susceptibilities
 block_gf_2t_t g_l, g_g;                     // Lesser and greater GF's to be calculated
 block_gf_2t_t g_ret, g_adv;                 // Retarded and advanced GF's to be calculated
 init_state const * initial_state = nullptr; // Initial state at t=t_min
 triqs::mpi::communicator comm;              // MPI communicator
 compute_gf_parameters_t compute_gf_params;  // Parameters of the last call to solve
 gf_mesh<retime> t_mesh;                     // 1D time mesh to use in calculations

 // Make g_ret and g_adv out of g_l and g_g
 void make_gf_ret_adv();

public:

 solver(gf_struct_t const& gf_struct, gf_struct_t const& chi_struct,
        std::pair<double,double> time_window, int n_t);

 /// Compute the Green's functions for given initial state and Hamiltonian
 TRIQS_WRAP_ARG_AS_DICT
 void compute_gf(compute_gf_parameters_t const& p);

 /// Get the initial state at t=t_min
 init_state const& get_initial_state() const {
  if(initial_state == nullptr) TRIQS_RUNTIME_ERROR << "Initial state is not set!";
  return *initial_state;
 }

 /// Set the initial state
 void set_initial_state(init_state const& initial_state) { this->initial_state = &initial_state; }

 /// Set of parameters used in the last call to compute_gf()
 compute_gf_parameters_t get_last_compute_gf_parameters() const { return compute_gf_params; }

 /// Lesser GF in real time
 block_gf_2t_t get_g_l() { return g_l; }

 /// Greater GF in real time
 block_gf_2t_t get_g_g() { return g_g; }

 /// Retarded GF in real time
 block_gf_2t_t get_g_ret() { return g_ret; }

 /// Advanced GF in real time
 block_gf_2t_t get_g_adv() { return g_adv; }
};

}
