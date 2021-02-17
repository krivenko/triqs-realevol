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

#include "wl_worker.hpp"

#include "array_utility.hpp"
#include "time_expr.hpp"
#include "time_interp.hpp"

namespace realevol {

//
// Constructor
//

template<std::size_t NPoints, typename HamiltonianType, typename TPointSelector>
wl_worker<NPoints, HamiltonianType, TPointSelector>::wl_worker(
  init_state const& initial_state,
  HamiltonianType const& h,
  double hbar,
  hilbert_space_structure<HamiltonianType> const& hss,
  std::vector<bool> const& is_static_sp,
  TPointSelector const& t_selector,
  gf_mesh<retime> const& t_mesh,
  h_interpolation HInterpol,
  lanczos_params_t const& lanczos_params
) :
  initial_state(initial_state),
  hss(hss),
  is_static_sp(is_static_sp),
  t_selector(t_selector),
  h(h, initial_state.get_fops(), initial_state.get_full_hs()),
  hbar(hbar),
  t_mesh(t_mesh),
  HInterpol(HInterpol),
  lanczos_params(lanczos_params) {
}

//
// Add a value to a GF element accessed by a std::array of indices.
// Too bad, TRIQS does not seemingly support this operation out of the box.
//

template<std::size_t... Ints>
void add_to_element_impl(time_container_view_t<sizeof...(Ints)> gf,
                         std::array<int, sizeof...(Ints)> const& t_indices,
                         dcomplex val,
                         std::index_sequence<Ints...>
                        ) {
  gf[{t_indices[Ints]...}] += val;
}

template<std::size_t NPoints>
void add_to_element(time_container_view_t<NPoints> gf,
                    std::array<int, NPoints> const& t_indices,
                    dcomplex val
                   ) {
  add_to_element_impl(gf, t_indices, val, std::make_index_sequence<NPoints>());
}

//
// Call time_point_selector<NPoints> with a mesh object and a std::array
// of time point indices.
//

template<std::size_t... Ints>
bool call_t_selector_impl(time_point_selector<sizeof...(Ints)> const& t_selector,
                          gf_mesh<retime> const& t_mesh,
                          std::array<int, sizeof...(Ints)> const& t_indices,
                          std::index_sequence<Ints...>
                         ) {
  return t_selector({t_mesh[t_indices[Ints]]...});
}

template<std::size_t NPoints>
bool call_t_selector(time_point_selector<NPoints> const& t_selector,
                     gf_mesh<retime> const& t_mesh,
                     std::array<int, NPoints> const& t_indices
                    ) {
  return call_t_selector_impl(t_selector,
                              t_mesh,
                              t_indices,
                              std::make_index_sequence<NPoints>());
}

//
// operator() implements a rather involved recursive algorithm.
// It evaluates a matrix element of the following form,
//
// <\psi_f| U(0, t_0) M_{N-1} U(t_0, t_1) ...
//                ... M_1 U(t_{N-2}, t_{N-1}) M_0 U(t_{N-1}, 0) |\psi_i>
//
// Here, |\psi_i> and |\psi_f> are initial and final states. They are elements of
// hss.sub_hilbert_spaces[0] and hss.sub_hilbert_spaces[N] respectively.
// Real time propagators U(t_{N-1}, 0), U(t_{N-2}, t_{N-1}), ..., U(0, t_0) act
// on states in subspaces hss.sub_hilbert_spaces[n] with n=0, ..., N.
// Monomial operators M_0, M_1, ..., M_{N-1} generate one-to-one connections
// between pairs of said subspaces.
//
// Let us define an auxiliary state <\phi(t_0)|
//
//   <\phi(t_0)| = <\psi_f| U(0, t_0),
//
// or, in the ket-vector notation,
//
//   |\phi(t_0)> = U(t_0, 0) |\psi_f>.
//
// We also define N+1 auxiliary vectors |\psi_n>,
//
//  |\psi_0> = U(t_{N-1}, 0) |\psi_i>
//  |\psi_1> = U(t_{N-2}, t_{N-1}) M_0 |\psi_0>
//  |\psi_2> = U(t_{N-3}, t_{N-2}) M_1 |\psi_1>
//
//  ...
//
//  |\psi_{N-1}> = U(t_0, t_1) M_{N-2} |\psi_{N-2}>
//  |\psi_N> = M_{N-1} |\psi_{N-1}>
//
// For a fixed sequence of time points t_0, t_1, ..., t_{N-1}, the matrix element
// in question equals <\phi(t_0)|\psi_N>
//
// The actual calculation proceeds as follows. In the outer loop over t_0,
// |\phi(t_0)> is incrementally calculated from |\psi_f>. At each iteration of this
// loop, |\psi_0> is set to |\psi_i> and do_inner_loop<0>() is called. The role
// of do_inner_loop<P>() is to calculate |\psi_{P+1}> for all values of the
// time point t_{N-1-P} and to recursively call do_inner_loop<P+1>() until
// |\psi_N> is reached.
//
// For instance, do_inner_loop<0>() runs a loop over t_{N-1} and at each iteration
// does the following:
//
// - Incrementally updates |\psi_0> -> U(t_{N-1}, 0) |\psi_0>
// - Computes |\psi_1> = M_0 |\psi_0>
// - Calls do_inner_loop<1>()
//
// Similarly, do_inner_loop<1>() iterates over t_{N-2} and
//
// - Incrementally updates |\psi_1> -> U(t_{N-2}, t_{N-1}) |\psi_1>
// - Computes |\psi_2> = M_1 |\psi_1>
// - Calls do_inner_loop<2>()
//
// This recursion terminates in do_inner_loop<N-1>() where
//
// - |\psi_{N-1}> is updated, |\psi_{N-1}> -> U(t_0, t_1) |\psi_{N-1}>
// - |\psi_N> is computed, |\psi_N> = M_{N-1} |\psi_{N-1}>
// - The final result is obtained as <\phi(t_0)|\psi_N>.
//

template<std::size_t NPoints /* N */,
         typename HamiltonianType,
         typename TPointSelector>
template<std::size_t Point>
void wl_worker<NPoints, HamiltonianType, TPointSelector>::do_inner_loop(
  std::array<propagator<HScalarType>, NPoints+1> const& props,
  std::array<op_with_map_t, NPoints> const& ops,
  std::array<state_on_subspace_t, NPoints+1> & psi,
  state_on_subspace_t & phi,
  dcomplex coeff,
  time_container_view_t<NPoints> result
) const {
  if constexpr(Point < NPoints - 1) { // Do a step of recursion

    int & t_index = t_indices[NPoints-1-Point];
    int t_index_prev = 0;

    // TODO: Get rid of this rewind operation if possible
    if constexpr(Point > 0)
      props[Point](psi[Point], t_indices[NPoints - Point], 0);

    for(t_index = 0; t_index < t_mesh.size(); ++t_index) {
      props[Point](psi[Point], t_index_prev, t_index);
      ops[Point].apply(psi[Point], psi[Point + 1]);

      if(!is_zero(dot_product(psi[Point + 1], psi[Point + 1])))
        do_inner_loop<Point + 1>(props, ops, psi, phi, coeff, result);

      t_index_prev = t_index;
    }
  } else {
    // The innermost loop is reached, terminating the recursion and
    // computing the scalar product.
    if(!call_t_selector(t_selector, t_mesh, t_indices)) return;

    if constexpr(NPoints == 1)
      props[NPoints - 1](psi[NPoints - 1], 0, t_indices[0]);
    else
      props[NPoints - 1](psi[NPoints - 1], t_indices[1], t_indices[0]);
    ops[NPoints - 1].apply(psi[NPoints - 1], psi[NPoints]);

    add_to_element(result, t_indices, coeff * dot_product(phi, psi[NPoints]));
  }
}

//
// Call operator
//

template<std::size_t NPoints, typename HamiltonianType, typename TPointSelector>
void wl_worker<NPoints, HamiltonianType, TPointSelector>::operator()(
  worldline_desc_t<NPoints> const& wl,
  time_container_view_t<NPoints> result
) const {

  // Propagators U(t_{N-1}, 0), U(t_{N-2}, t_{N-1}), ..., U(0, t_0)
  auto props = map_array<propagator<HScalarType>>([&](long sp) {
    return propagator<HScalarType>(
      h,
      hss.sub_hilbert_spaces[sp],
      t_mesh,
      hbar,
      is_static_sp[sp],
      HInterpol,
      lanczos_params.min_matrix_size,
      lanczos_params.get_energy_tol(sp),
      lanczos_params.get_max_krylov_dim(sp)
    );
  }, wl.sp_indices);

  auto const& fops = initial_state.get_fops();
  auto const& full_hs = initial_state.get_full_hs();

  // Monomial operators M_0, M_1, ...
  auto ops = map_array<op_with_map_t>([&](auto const& m) {
    return op_with_map_t(
      many_body_operator(1, m),
      fops,
      full_hs,
      std::move(hss.make_monomial_connections(m)),
      &hss.sub_hilbert_spaces
    );
  }, wl.M);

  auto const& wst = initial_state.get_weighted_states()[wl.weighted_state_index];
  auto const& subspaces = hss.sub_hilbert_spaces;

  // Overall prefactor of this world line
  dcomplex coeff = wl.factor * wst.weight;

  // Initial state |\psi_i>
  auto psi_i = project<state_on_subspace_t>(wst.state, subspaces[wl.sp_indices[0]]);

  // State |\phi(t = 0)>
  auto phi = project<state_on_subspace_t>(wst.state, subspaces[wl.sp_indices[NPoints]]);

  // Intermediate states |\psi_n>
  auto psi = map_array<state_on_subspace_t>([&](long sp) {
    return state_on_subspace_t(hss.sub_hilbert_spaces[sp]);
  }, wl.sp_indices);

  // Outer loop over t_0
  int t_index_prev = 0;
  int & t_index = t_indices[0];
  for(t_index = 0; t_index < t_mesh.size(); ++t_index) {

    // Update |\phi(t)>
    props[NPoints](phi, t_index_prev, t_index);

    // Enter the first inner loop
    psi[0] = psi_i;
    do_inner_loop<0>(props, ops, psi, phi, coeff, result);

    t_index_prev = t_index;
  }
}

//
// Explicit instantiations
//

using time_expr_operator_t = realevol::operators::many_body_operator_generic<time_expr>;
using time_interp_operator_t = realevol::operators::many_body_operator_generic<time_interp>;

// Compute expectation values
template class wl_worker<1, time_expr_operator_t, time_point_selector<1>>;
template class wl_worker<1, time_interp_operator_t, time_point_selector<1>>;

// Compute 2-point correlators
template class wl_worker<2, time_expr_operator_t, time_point_selector<2>>;
template class wl_worker<2, time_interp_operator_t, time_point_selector<2>>;

// Compute GFs and susceptibilities
template class wl_worker<2, time_expr_operator_t, time_point_selector_lower_triangle>;
template class wl_worker<2, time_interp_operator_t, time_point_selector_lower_triangle>;

// Compute 3-point correlators
template class wl_worker<3, time_expr_operator_t, time_point_selector<3>>;
template class wl_worker<3, time_interp_operator_t, time_point_selector<3>>;

} // namespace realevol
