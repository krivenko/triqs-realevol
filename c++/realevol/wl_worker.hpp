/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

#include "types.hpp"
#include "init_state.hpp"
#include "worldlines.hpp"
#include "propagator.hpp"
#include "time_point_selector.hpp"
#include "dynamical_trace.hpp"

namespace realevol {

// Evaluate worldlines
template<typename HamiltonianType>
class wl_worker {

 using HScalarType = typename HamiltonianType::scalar_t;

 init_state const& initial_state;
 hilbert_space_structure<HamiltonianType> const& hss;
 HamiltonianType const& h_;
 op_on_subspace_t<HScalarType> h;
 double hbar;

 std::vector<sub_hilbert_space> const& subspaces;
 std::vector<bool> const& is_static_sp;
 mutable time_point_selector<2> t_selector;
 std::vector<std::vector<int>> c_conn, cdag_conn, n_conn;

 int lanczos_min_matrix_size;
 std::map<long,double> const& lanczos_gs_energy_tol;
 std::map<long,int> const& lanczos_max_krylov_dim;

public:

 wl_worker(init_state const& initial_state, HamiltonianType const& h, double hbar,
           hilbert_space_structure<HamiltonianType> const& hss,
           std::vector<bool> const& is_static_sp,
           time_point_selector<2> t_selector,
           int lanczos_min_matrix_size,
           std::map<long,double> const& lanczos_gs_energy_tol,
           std::map<long,int> const& lanczos_max_krylov_dim
          ) :
 initial_state(initial_state), hss(hss),
 h_(h),
 h(h, initial_state.get_fops(), initial_state.get_full_hs()), hbar(hbar),
 subspaces(hss.sub_hilbert_spaces),
 is_static_sp(is_static_sp),
 t_selector(std::move(t_selector)),
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
 void operator()(worldline_desc_t<2> const& wl, gf_2t_t & obs,
                 std::integral_constant<h_interpolation,HInterpol>) const {

  auto t_mesh = std::get<0>(obs.mesh());

  auto lanczos_params = lanczos_params_t{lanczos_min_matrix_size,
                                         lanczos_gs_energy_tol,
                                         lanczos_max_krylov_dim};

  t_selector.set_swap_t_tp(wl.observable == worldline_desc_t<2>::LesserGf);
  dynamical_trace<2, HamiltonianType> trace(
    initial_state,
    h_,
    hbar,
    hss,
    is_static_sp,
    t_selector,
    t_mesh,
    HInterpol,
    lanczos_params
  );

  auto result = correlator_2t_container_t{{t_mesh, t_mesh}};
  trace(wl, result);

  for(auto t1 : t_mesh) {
    for(auto t2 : t_mesh) {
      (wl.observable == worldline_desc_t<2>::LesserGf ?
        obs[{t2, t1}] : obs[{t1, t2}])(wl.inner_index1, wl.inner_index2)
          += result[{t1, t2}];
    }
  }
 }

};

}
