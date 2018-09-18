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
 std::map<long,double> const& lanczos_gs_energy_tol;
 std::map<long,int> const& lanczos_max_krylov_dim;

 std::vector<sub_hilbert_space> const& subspaces;
 std::vector<bool> const& is_static_sp;
 std::vector<std::vector<int>> c_conn, cdag_conn, n_conn;

public:

 wl_worker(init_state const& initial_state, operator_t const& h, double hbar,
           hilbert_space_structure const& hss,
           std::vector<bool> const& is_static_sp,
           int lanczos_min_matrix_size,
           std::map<long,double> const& lanczos_gs_energy_tol,
           std::map<long,int> const& lanczos_max_krylov_dim
          ) :
 initial_state(initial_state), hss(hss), subspaces(hss.sub_hilbert_spaces),
 h(h, initial_state.get_fops(), initial_state.get_full_hs()), hbar(hbar),
 is_static_sp(is_static_sp),
 lanczos_min_matrix_size(lanczos_min_matrix_size),
 lanczos_gs_energy_tol(lanczos_gs_energy_tol),
 lanczos_max_krylov_dim(lanczos_max_krylov_dim) {

  auto const& fops = initial_state.get_fops();
  c_conn.resize(fops.size(Fermion), std::vector<int>(subspaces.size()));
  cdag_conn.resize(fops.size(Fermion), std::vector<int>(subspaces.size()));
  n_conn.resize(fops.size(Fermion), std::vector<int>(subspaces.size()));

  for(int linear_index : range(fops.size(Fermion))) {
   for(int spn : range(subspaces.size())) {
    int c_conn_sp = hss.annihilation_connection[Fermion](linear_index, spn);
    c_conn[linear_index][spn] = c_conn_sp;
    cdag_conn[linear_index][spn] = hss.creation_connection[Fermion](linear_index, spn);
    n_conn[linear_index][spn] = (c_conn_sp == -1 ? -1 :
     hss.creation_connection[Fermion](linear_index, c_conn_sp));
   }
  }
 }

 // Do the actual observable calculation for one worldline
 template<h_interpolation HInterpol>
 void operator()(worldline_desc_t const& wl, gf_2t_t & obs,
                 std::integral_constant<h_interpolation,HInterpol>) const {

  auto const& left_hs = subspaces[wl.left_sp_index];
  auto const& middle_hs = subspaces[wl.middle_sp_index];
  auto const& right_hs = subspaces[wl.right_sp_index];

  auto get_gs_energy_tol = [this](int spn) {
   auto it = lanczos_gs_energy_tol.find(spn);
   return it != lanczos_gs_energy_tol.end() ? it->second : -1;
  };
  auto get_max_krylov_dim = [this](int spn) {
   auto it = lanczos_max_krylov_dim.find(spn);
   return it != lanczos_max_krylov_dim.end() ? it->second : -1;
  };

  propagator left_prop(h, left_hs, hbar, is_static_sp[wl.left_sp_index], HInterpol,
                       lanczos_min_matrix_size,
                       get_gs_energy_tol(wl.left_sp_index),
                       get_max_krylov_dim(wl.left_sp_index)
                      );
  propagator middle_prop(h, middle_hs, hbar, is_static_sp[wl.middle_sp_index], HInterpol,
                         lanczos_min_matrix_size,
                         get_gs_energy_tol(wl.middle_sp_index),
                         get_max_krylov_dim(wl.middle_sp_index)
                        );
  propagator right_prop(h, right_hs, hbar, is_static_sp[wl.right_sp_index], HInterpol,
                        lanczos_min_matrix_size,
                        get_gs_energy_tol(wl.right_sp_index),
                        get_max_krylov_dim(wl.right_sp_index)
                       );

  auto const& fops = initial_state.get_fops();
  auto const& full_hs = initial_state.get_full_hs();

  auto const& wst = initial_state.get_weighted_states()[wl.weighted_state_index];

  using op_with_map_t = imperative_operator<sub_hilbert_space,dcomplex,true>;
  op_with_map_t A; // operator connecting middle_hs to left_hs
  op_with_map_t B; // operator connecting right_hs to middle_hs
  // 1d time meshes associated with operators A and B
  gf_mesh<retime> A_mesh, B_mesh;
  dcomplex coeff = wst.weight;

  switch(wl.observable) {
   case worldline_desc_t::GreaterGf:
    A = op_with_map_t(c(wl.index1), fops, full_hs, c_conn[fops[wl.index1]], &subspaces);
    B = op_with_map_t(c_dag(wl.index2), fops, full_hs, cdag_conn[fops[wl.index2]], &subspaces);
    A_mesh = std::get<0>(obs.mesh());
    B_mesh = std::get<1>(obs.mesh());
    coeff *= -1_j / hbar;
    break;
   case worldline_desc_t::LesserGf:
    A = op_with_map_t(c_dag(wl.index2), fops, full_hs, cdag_conn[fops[wl.index2]], &subspaces);
    B = op_with_map_t(c(wl.index1), fops, full_hs, c_conn[fops[wl.index1]], &subspaces);
    A_mesh = std::get<1>(obs.mesh());
    B_mesh = std::get<0>(obs.mesh());
    coeff *= 1_j / hbar;
    break;
   case worldline_desc_t::Susceptibility:
    A = op_with_map_t(n(wl.index1), fops, full_hs, n_conn[fops[wl.index1]], &subspaces);
    B = op_with_map_t(n(wl.index2), fops, full_hs, n_conn[fops[wl.index2]], &subspaces);
    A_mesh = std::get<0>(obs.mesh());
    B_mesh = std::get<1>(obs.mesh());
    coeff *= -1_j / hbar;
    break;
  }

  // Worldline contribution is evaluated as
  // coeff * <bra_st(*A_it)| A U(*A_it,*B_it) B |ket_st(*B_it)>

  // <bra_st(*A_it)| = <psi_0| U(0,*A_it)
  // |bra_st(*A_it)> = U(*A_it,0) |psi_0>
  auto bra_st = project<state_on_subspace_t>(wst.state, left_hs);
  auto A_it_prev = A_mesh.begin();
  for(auto A_it = A_mesh.begin(); A_it != A_mesh.end(); A_it_prev = A_it++) {
   left_prop(bra_st, A_it_prev, A_it);

   // |ket_st(*B_it)> = U(*B_it,0)|psi_0>
   auto ket_st = project<state_on_subspace_t>(wst.state, right_hs);
   auto B_it_prev = B_mesh.begin();
   for(auto B_it = B_mesh.begin(); B_it != B_mesh.end(); B_it_prev = B_it++) {
    right_prop(ket_st, B_it_prev, B_it);
 
    // |middle_st(*A_it)> = U(*A_it,*B_it) B |ket_st(*B_it)>
    auto middle_st = state_on_subspace_t(middle_hs);
    B.apply(ket_st, middle_st);
    
    // |a_st(*A_it)> = A |middle_st(*A_it)>
    auto a_st = state_on_subspace_t(left_hs);
    middle_prop(middle_st, B_it, A_it);
    A.apply(middle_st, a_st);

    (wl.observable == worldline_desc_t::LesserGf ?
     obs[{*B_it,*A_it}] : obs[{*A_it,*B_it}])
     (wl.inner_index1, wl.inner_index2) += coeff * dot_product(bra_st, a_st);
   }
  }
 }

};

}
