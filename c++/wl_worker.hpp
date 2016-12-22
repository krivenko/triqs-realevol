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

#include <vector>
#include <type_traits>

#include <triqs/utility/signal_handler.hpp>

#include "common.hpp"
#include "init_state.hpp"
#include "worldlines.hpp"
#include "propagator.hpp"

namespace realevol {

// Evaluate worldlines
class wl_worker {

 init_state const& initial_state;
 hilbert_space_structure const& hss;
 op_on_subspace_t h;
 double hbar;

 int lanczos_min_matrix_size;

 std::vector<sub_hilbert_space> const& subspaces;
 std::vector<std::vector<int>> c_conn, cdag_conn, n_conn;

public:

 wl_worker(init_state const& initial_state, operator_t const& h, double hbar,
           hilbert_space_structure const& hss,
           int lanczos_min_matrix_size) :
 initial_state(initial_state), hss(hss), subspaces(hss.sub_hilbert_spaces),
 h(h, initial_state.get_fops(), initial_state.get_full_hs()), hbar(hbar),
 lanczos_min_matrix_size(lanczos_min_matrix_size) {

  auto const& fops = initial_state.get_fops();
  c_conn.resize(fops.size(Fermion), std::vector<int>(subspaces.size()));
  cdag_conn.resize(fops.size(Fermion), std::vector<int>(subspaces.size()));
  n_conn.resize(fops.size(Fermion), std::vector<int>(subspaces.size()));
  for(int linear_index : range(fops.size(Fermion))) {
   for(int spn : range(subspaces.size())) {
    c_conn[linear_index][spn]
     = hss.annihilation_connection[Fermion](linear_index, spn);
    cdag_conn[linear_index][spn]
     = hss.creation_connection[Fermion](linear_index, spn);
    n_conn[linear_index][spn] = cdag_conn[linear_index][c_conn[linear_index][spn]];
   }
  }
 }

 // Do the actual observable calculation for one worldline
 template<h_interpolation HInterpol>
 void operator()(worldline_desc_t const& wl, block_gf_2t_t & obs,
                 std::integral_constant<h_interpolation,HInterpol>) const {

  auto const& left_hs = subspaces[wl.left_sp_index];
  auto const& middle_hs = subspaces[wl.middle_sp_index];
  auto const& right_hs = subspaces[wl.right_sp_index];

  propagator<HInterpol> left_prop(h, left_hs, hbar, lanczos_min_matrix_size);
  propagator<HInterpol> middle_prop(h, middle_hs, hbar, lanczos_min_matrix_size);
  propagator<HInterpol> right_prop(h, right_hs, hbar, lanczos_min_matrix_size);

  auto const& fops = initial_state.get_fops();
  auto const& full_hs = initial_state.get_full_hs();

  auto const& wst = initial_state.get_weighted_states()[wl.weighted_state_index];

  using op_with_map_t = imperative_operator<sub_hilbert_space,dcomplex,true>;
  op_with_map_t left_op, right_op;
  dcomplex coeff = wst.weight;

  if(wl.observable == worldline_desc_t::GreensFunction) { // GF
   if(wl.is_greater) {
    left_op = op_with_map_t(c(wl.index1), fops, full_hs,
                            c_conn[fops[wl.index1]], &subspaces);
    right_op = op_with_map_t(c_dag(wl.index2), fops, full_hs,
                             cdag_conn[fops[wl.index2]], &subspaces);
    coeff *= -1_j / hbar;
   } else {
    left_op = op_with_map_t(c_dag(wl.index2), fops, full_hs,
                            cdag_conn[fops[wl.index2]], &subspaces);
    right_op = op_with_map_t(c(wl.index1), fops, full_hs,
                             c_conn[fops[wl.index1]], &subspaces);
    coeff *= 1_j / hbar;
   }
  } else { // Susceptibility
   if(wl.is_greater) {
    left_op = op_with_map_t(n(wl.index1), fops, full_hs,
                            n_conn[fops[wl.index1]], &subspaces);
    right_op = op_with_map_t(n(wl.index2), fops, full_hs,
                             n_conn[fops[wl.index2]], &subspaces);
   } else {
    left_op = op_with_map_t(n(wl.index2), fops, full_hs,
                            n_conn[fops[wl.index2]], &subspaces);
    right_op = op_with_map_t(n(wl.index1), fops, full_hs,
                             n_conn[fops[wl.index1]], &subspaces);
   }
   coeff *= -1_j / hbar;
  }

  auto bra_st = project<state_on_subspace_t>(wst.state, left_hs);
  auto right_st = project<state_on_subspace_t>(wst.state, right_hs);
  auto middle_st = state_on_subspace_t(middle_hs);
  auto left_st = state_on_subspace_t(middle_hs);

  gf_2t_view obs_block = obs[wl.block_index];

  gf_mesh<retime> left_t_mesh = (wl.is_greater ? std::get<0>(obs_block.mesh()) :
                                                 std::get<1>(obs_block.mesh()));
  gf_mesh<retime> right_t_mesh = (wl.is_greater ? std::get<1>(obs_block.mesh()) :
                                                  std::get<0>(obs_block.mesh()));

  auto right_it_prev = right_t_mesh.begin();
  for(auto right_it = right_t_mesh.begin(); right_it != right_t_mesh.end();
      right_it_prev = right_it++) {
   right_prop(right_st, right_it_prev, right_it);
   middle_st = right_op(right_st);
   auto left_it_prev = right_it;
   for(auto left_it = left_t_mesh.begin(); left_it != left_t_mesh.end();
       left_it_prev = left_it++) {
    middle_prop(middle_st, left_it_prev, left_it);
    left_st = left_op(middle_st);
    left_prop(left_st, left_it, left_t_mesh.begin());

    (wl.is_greater ? obs_block[{*left_it,*right_it}] : obs_block[{*right_it,*left_it}])
     (wl.inner_index1, wl.inner_index2) += coeff * dot_product(bra_st, left_st);
   }
  }
 }

};

}
