/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <string>
#include <vector>
#include <map>
#include <utility>

#include <mpi/mpi.hpp>

#include "triqs/operators/many_body_operator.hpp"
#include "compute_2t_obs_parameters.hpp"
#include "solver_parameters.hpp"


namespace realevol {

using namespace triqs::gfs;

class solver {

 gf_struct_t gf_struct;                              // Block structure of the Green functions
 chi_indices_t chi_indices;                          // Indices of susceptibilities
 block_gf_2t_t g_l, g_g;                             // Lesser and greater GF's to be calculated
 gf_2t_t chi;                                        // Susceptibility to be calculated
 init_state const * initial_state = nullptr;         // Initial state at t=0
 mpi::communicator comm;                             // MPI communicator
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
 void set_initial_state(init_state const& state) { this->initial_state = &state; }

 /// Set of parameters used in the last call to compute_2t_obs()
 compute_2t_obs_parameters_t get_last_compute_2t_obs_parameters() const {
   return compute_2t_obs_params;
 }

 /// Lesser GF in real time
 block_gf_2t_view get_g_l() { return g_l; }

 /// Greater GF in real time
 block_gf_2t_view get_g_g() { return g_g; }

 /// Susceptibility in real time
 gf_2t_view get_chi() { return chi; }

};

//
// Expectation values
//

// Compute expectation value of operator 'op' as a function of time
template<typename HamiltonianType>
expectval_container_t compute_expectval(static_operator_t const& op,
                                        init_state const& initial_state,
                                        HamiltonianType const& h,
                                        mesh_t_t const& t_mesh,
                                        solver_parameters_t const& params,
                                        mpi::communicator const& comm = {});

//
// 2-point correlator
//

// Compute a 2-point correlator of operators 'op1' and 'op2' with their time
// arguments defined on a Cartesian product 't_mesh' x 't_mesh'.
template<typename HamiltonianType>
correlator_2t_container_t compute_correlator_2t(static_operator_t const& op1,
                                                static_operator_t const& op2,
                                                init_state const& initial_state,
                                                HamiltonianType const& h,
                                                mesh_t_t const& t_mesh,
                                                solver_parameters_t const& params,
                                                mpi::communicator const& comm = {});

//
// 3-point correlator
//

// Compute a 3-point correlator of operators 'op1', 'op2' and 'op3' with their time
// arguments defined on a Cartesian product 't_mesh' x 't_mesh' x 't_mesh'.
template<typename HamiltonianType>
correlator_3t_container_t compute_correlator_3t(static_operator_t const& op1,
                                                static_operator_t const& op2,
                                                static_operator_t const& op3,
                                                init_state const& initial_state,
                                                HamiltonianType const& h,
                                                mesh_t_t const& t_mesh,
                                                solver_parameters_t const& params,
                                                mpi::communicator const& comm = {});

//
// Lesser GF
//

// Compute the lesser GF defined on a Cartesian product 't_mesh' x 't_mesh'.
template<typename HamiltonianType>
block_gf_2t_t compute_g_l(gf_struct_t const& gf_struct,
                          init_state const& initial_state,
                          HamiltonianType const& h,
                          mesh_t_t const& t_mesh,
                          solver_parameters_t const& params,
                          mpi::communicator const& comm = {});

//
// Greater GF
//

// Compute the greater GF defined on a Cartesian product 't_mesh' x 't_mesh'.
template<typename HamiltonianType>
block_gf_2t_t compute_g_g(gf_struct_t const& gf_struct,
                          init_state const& initial_state,
                          HamiltonianType const& h,
                          mesh_t_t const& t_mesh,
                          solver_parameters_t const& params,
                          mpi::communicator const& comm = {});

//
// Susceptibility
//

// Compute the susceptibility defined on a Cartesian product 't_mesh' x 't_mesh'.
template<typename HamiltonianType>
block_gf_2t_t compute_chi(gf_struct_t const& gf_struct,
                          init_state const& initial_state,
                          HamiltonianType const& h,
                          mesh_t_t const& t_mesh,
                          solver_parameters_t const& params,
                          mpi::communicator const& comm = {});
}
