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

/////////////////
/// diag_1d_t ///
/////////////////

inline diag_1d_t::diag_1d_t(propagator const& prop, sub_hilbert_space const& sp) :
 prop(prop), tmp_st(sp) {
 tmp_st(0) = 1;
}

inline void diag_1d_t::diag(double t, double dt) {
 eigenvalue = prop.apply_h(tmp_st, t, dt)(0).real();
}

/////////////////////
/// diag_lapack_t ///
/////////////////////

inline diag_lapack_t::diag_lapack_t(propagator const& prop, sub_hilbert_space const& sp) :
 prop(prop), from_st(sp), to_st(sp), workspace(prop.N, prop.N) {
}

inline void diag_lapack_t::diag(double t, double dt) {
 for(long i : range(prop.N)) {
  from_st(i) = 1.0;
  to_st = prop.apply_h(from_st, t, dt);
  workspace(range(),i) = to_st.amplitudes();
  from_st(i) = 0;
 }
 eig = linalg::eigenelements_in_place(&workspace);
}

//////////////////
/// propagator ///
//////////////////s

propagator::propagator(op_on_subspace_t const& h, sub_hilbert_space const& sp,
                       double hbar, bool is_static_h, h_interpolation h_interpol,
                       long lanczos_min_matrix_size) :
 h(h), h_coeff(-1_j / hbar), N(sp.size()), is_static_h(is_static_h), h_interpol(h_interpol),
 ed_solver(N == 1 ? OneD : (N < lanczos_min_matrix_size ? LAPACK : Lanczos)) {
 switch(ed_solver) {
  case OneD:
  diag_1d = std::make_unique<diag_1d_t>(*this, sp);
  if(is_static_h) diag_1d->diag(0,0);
  break;
  case LAPACK:
  diag_lapack = std::make_unique<diag_lapack_t>(*this, sp);
  if(is_static_h) diag_lapack->diag(0,0);
  break;
  case Lanczos:
  // TODO
  break;
 }
}

void propagator::operator()(state_on_subspace_t & st,
                            time_it_t const& t_start,
                            time_it_t const& t_end) const {
 if(t_end == t_start) return; // No propagation needed

 bool forward = *t_end > *t_start;
 dcomplex c = (forward ? 1 : -1) * h_coeff;

 if(is_static_h) {
  if(forward) propagate(st, t_start, t_end, c); // Forward propagation (static)
  else        propagate(st, t_end, t_start, c); // Backward propagation (static)
 } else {
  if(forward) { // Forward propagation
   time_it_t t_next = t_start; ++t_next;
   for(auto t = t_start; t != t_end; ++t, ++t_next) propagate(st, t, t_next, c);
  } else { // Backward propagation
   time_it_t t_next = t_end; ++t_next;
   for(auto t = t_end; t != t_start; ++t, ++t_next) propagate(st, t, t_next, c);
  }
 }
}

inline void propagator::propagate(state_on_subspace_t & st,
                                  time_it_t const& t_min,
                                  time_it_t const& t_max,
                                  dcomplex c) const {
 switch(ed_solver) {
  case OneD:    return propagate(st, t_min, t_max, c, ed_type<OneD>());
  case LAPACK:  return propagate(st, t_min, t_max, c, ed_type<LAPACK>());
  case Lanczos: return propagate(st, t_min, t_max, c, ed_type<Lanczos>());
 }
}

/// 1D propagation
inline void propagator::propagate(state_on_subspace_t & st,
                                  time_it_t const& t_min, time_it_t const& t_max,
                                  dcomplex c, ed_type<OneD>) const {
 double dt = *t_max - *t_min;
 if(!is_static_h) diag_1d->diag(*t_min, dt);
 st *= std::exp(c * dt * diag_1d->eigenvalue);
}

/// LAPACK propagation
inline void propagator::propagate(state_on_subspace_t & st,
                                  time_it_t const& t_min, time_it_t const& t_max,
                                  dcomplex c, ed_type<LAPACK>) const {
 double dt = *t_max - *t_min;
 if(!is_static_h) diag_lapack->diag(*t_min, dt);
 auto & psi = st.amplitudes();
 psi = conj(diag_lapack->eig.second) * psi;
 for(long i : range(N)) psi(i) *= std::exp(c * dt * diag_lapack->eig.first(i));
 psi = diag_lapack->eig.second.transpose() * psi;
}

inline void propagator::propagate(state_on_subspace_t & st,
                                  time_it_t const& t_min, time_it_t const& t_max,
                                  dcomplex c, ed_type<Lanczos>) const {
 double dt = *t_max - *t_min;
 // TODO
}

inline state_on_subspace_t propagator::apply_h(state_on_subspace_t const& st, double t, double dt) const {
 switch(h_interpol) {
  case Rectangle: return h(st, t+dt/2);
  case Trapezoid: return (h(st, t) + h(st, t+dt)) / 2.0;
  case Simpson:   return (h(st, t) + 4*h(st, t+dt/2) + h(st, t+dt)) / 6.0;
 }
}

}
