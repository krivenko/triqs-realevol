/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2023, I. Krivenko, M. Danilov, P. Kubiczek
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

#include "../time_expr.hpp"

#include "impl.hxx"

namespace realevol {

template
std::vector<expectval_container_t>
compute_expectval(std::vector<static_operator_t> const& ops,
                  init_state const& initial_state,
                  time_expr_operator_t const& h,
                  mesh_t_t const& t_mesh,
                  solver_parameters_t<1> const& params,
                  mpi::communicator const& comm);

template
expectval_container_t compute_expectval(static_operator_t const& op,
                                        init_state const& initial_state,
                                        time_expr_operator_t const& h,
                                        mesh_t_t const& t_mesh,
                                        solver_parameters_t<1> const& params,
                                        mpi::communicator const& comm);

template
std::vector<correlator_2t_container_t>
compute_correlator_2t(std::vector<std::array<static_operator_t, 2>> ops,
                      init_state const& initial_state,
                      time_expr_operator_t const& h,
                      mesh_t_t const& t_mesh,
                      solver_parameters_t<2> const& params,
                      mpi::communicator const& comm);

template
correlator_2t_container_t compute_correlator_2t(static_operator_t const& op1,
                                                static_operator_t const& op2,
                                                init_state const& initial_state,
                                                time_expr_operator_t const& h,
                                                mesh_t_t const& t_mesh,
                                                solver_parameters_t<2> const& params,
                                                mpi::communicator const& comm);

template
std::vector<correlator_3t_container_t>
compute_correlator_3t(std::vector<std::array<static_operator_t, 3>> ops,
                      init_state const& initial_state,
                      time_expr_operator_t const& h,
                      mesh_t_t const& t_mesh,
                      solver_parameters_t<3> const& params,
                      mpi::communicator const& comm);

template
correlator_3t_container_t compute_correlator_3t(static_operator_t const& op1,
                                                static_operator_t const& op2,
                                                static_operator_t const& op3,
                                                init_state const& initial_state,
                                                time_expr_operator_t const& h,
                                                mesh_t_t const& t_mesh,
                                                solver_parameters_t<3> const& params,
                                                mpi::communicator const& comm);

template
block_gf_2t_t compute_g_l(gf_struct_t const& gf_struct,
                          init_state const& initial_state,
                          time_expr_operator_t const& h,
                          mesh_t_t const& t_mesh,
                          solver_parameters_t<2> const& params,
                          mpi::communicator const& comm);

template
block_gf_2t_t compute_g_g(gf_struct_t const& gf_struct,
                          init_state const& initial_state,
                          time_expr_operator_t const& h,
                          mesh_t_t const& t_mesh,
                          solver_parameters_t<2> const& params,
                          mpi::communicator const& comm);

template
gf_2t_t compute_chi(chi_indices_t const& chi_indices,
                    init_state const& initial_state,
                    time_expr_operator_t const& h,
                    mesh_t_t const& t_mesh,
                    solver_parameters_t<2> const& params,
                    mpi::communicator const& comm);

} // namespace realevol
