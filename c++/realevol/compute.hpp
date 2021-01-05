/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2021, I. Krivenko, M. Danilov, P. Kubiczek
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

#include "types.hpp"
#include "init_state.hpp"
#include "solver_parameters.hpp"

namespace realevol {

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
gf_2t_t compute_chi(chi_indices_t const& chi_indices,
                    init_state const& initial_state,
                    HamiltonianType const& h,
                    mesh_t_t const& t_mesh,
                    solver_parameters_t const& params,
                    mpi::communicator const& comm = {});
}
