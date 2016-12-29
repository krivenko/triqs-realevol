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

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

#include <cmath>
#include <triqs/arrays/linalg/eigenelements.hpp>

#include "common.hpp"
#include "propagator.hpp"

namespace realevol {

using namespace triqs::arrays;

template<h_interpolation HInterpol>
propagator<HInterpol>::propagator(op_on_subspace_t const& h, sub_hilbert_space const& sp,
                               double hbar, int lanczos_min_matrix_size) :
 h(h), h_coeff(-1_j / hbar), N(sp.size()) {
 if(N == 1) {
  propagate = &propagator::propagate_1d;
  state_1d = state_on_subspace_t(sp);
  state_1d(0) = 1;
 } else if(N < lanczos_min_matrix_size) {
  propagate = &propagator::propagate_lapack;
  lapack_workspace = matrix<dcomplex>(N,N);
  lapack_st_from = state_on_subspace_t(sp);
  lapack_st_to = state_on_subspace_t(sp);
 } else {
  propagate = &propagator::propagate_lanczos;
 }
}

template<h_interpolation HInterpol>
void propagator<HInterpol>::operator()(state_on_subspace_t & st,
                                    time_it_t const& t_start, time_it_t const& t_end) const {
 if(t_end == t_start)       // No propagation needed
  return;
 else if(*t_end > *t_start) // Forward propagation
  (this->*propagate)(st, t_start, t_end, true);
 else                       // Backward propagation
  (this->*propagate)(st, t_end, t_start, false);
}

template<h_interpolation HInterpol>
void propagator<HInterpol>::propagate_1d(state_on_subspace_t & st,
                                      time_it_t t, time_it_t const& t_max, bool forward) const {
 dcomplex c = (forward ? 1 : -1) * h_coeff;
 for(; t != t_max; ++t) {
  auto t_next = t; ++t_next;
  double dt = *t_next - *t;
  st *= std::exp(c * dt * apply_h(state_1d, *t, dt)(0));
 }
}

template<h_interpolation HInterpol> void
propagator<HInterpol>::propagate_lapack(state_on_subspace_t & st,
                                     time_it_t t, time_it_t const& t_max, bool forward) const {

 dcomplex c = (forward ? 1 : -1) * h_coeff;
 auto const& hs = st.get_hilbert();

 for(; t != t_max; ++t) {
  auto t_next = t; ++t_next;
  double dt = *t_next - *t;

  for(long i : range(N)) {
   lapack_st_from(i) = 1.0;
   lapack_st_to = apply_h(lapack_st_from, *t, dt);
   lapack_workspace(range(),i) = lapack_st_to.amplitudes();
   lapack_st_from(i) = 0;
  }

  auto & psi = st.amplitudes();
  auto eig = linalg::eigenelements_in_place(&lapack_workspace);
  psi = conj(eig.second) * psi;
  for(long i : range(N)) psi(i) *= std::exp(c * dt * eig.first(i));
  psi = eig.second.transpose() * psi;
 }
}

template<h_interpolation HInterpol> void
propagator<HInterpol>::propagate_lanczos(state_on_subspace_t & st,
                                         time_it_t t, time_it_t const& t_max, bool forward) const {
 // TODO
}


template<> inline state_on_subspace_t
propagator<Rectangle>::apply_h(state_on_subspace_t const& st, double t, double dt) const {
 return h(st, t+dt/2);
}

template<> inline state_on_subspace_t
propagator<Trapezoid>::apply_h(state_on_subspace_t const& st, double t, double dt) const {
 return (h(st, t) + h(st, t+dt)) / 2.0;
}

template<> inline state_on_subspace_t
propagator<Simpson>::apply_h(state_on_subspace_t const& st, double t, double dt) const {
 return (h(st, t) + 4*h(st, t+dt/2) + h(st, t+dt)) / 6.0;
}

// Explicit instantiation
template class propagator<Rectangle>;
template class propagator<Trapezoid>;
template class propagator<Simpson>;

}
