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

////////////////
/// inst_h_t ///
////////////////

inline inst_h_t::inst_h_t(op_on_subspace_t const& h, bool is_static, h_interpolation h_interpol) :
 is_static(is_static), h(h), h_interpol(is_static ? Rectangle : h_interpol), t(0), dt(0) {
}

inline state_on_subspace_t inst_h_t::operator()(state_on_subspace_t const& st) const {
 switch(h_interpol) {
  case Rectangle: return h(st, t+dt/2);
  case Trapezoid: return (h(st, t) + h(st, t+dt)) / 2.0;
  case Simpson:   return (h(st, t) + 4*h(st, t+dt/2) + h(st, t+dt)) / 6.0;
 }
}

/////////////////
/// prop_1d_t ///
/////////////////

inline prop_1d_t::prop_1d_t(inst_h_t && inst_h_, sub_hilbert_space const& sp) :
 inst_h(std::move(inst_h_)), tmp_st(sp) {
 tmp_st(0) = 1;
 if(inst_h.is_static) diag(0,0);
}

inline void prop_1d_t::diag(double t, double dt) {
 inst_h.t = t;
 inst_h.dt = dt;
 eigenvalue = inst_h(tmp_st)(0).real();
}

inline void prop_1d_t::operator()(state_on_subspace_t & st,
                                  time_it_t const& t_min, time_it_t const& t_max,
                                  dcomplex c) {
 double dt = *t_max - *t_min;
 if(!inst_h.is_static) diag(*t_min, dt);
 st *= std::exp(c * dt * eigenvalue);
}

/////////////////////
/// prop_lapack_t ///
/////////////////////

inline prop_lapack_t::prop_lapack_t(inst_h_t && inst_h_, sub_hilbert_space const& sp) :
 inst_h(std::move(inst_h_)), from_st(sp), to_st(sp), workspace(sp.size(), sp.size()) {
 if(inst_h.is_static) diag(0,0);
}

inline void prop_lapack_t::diag(double t, double dt) {
 inst_h.t = t;
 inst_h.dt = dt;
 for(long i : range(from_st.size())) {
  from_st(i) = 1.0;
  to_st = inst_h(from_st);
  workspace(range(),i) = to_st.amplitudes();
  from_st(i) = 0;
 }
 eig = linalg::eigenelements_in_place(&workspace);
}

inline void prop_lapack_t::operator()(state_on_subspace_t & st,
                                      time_it_t const& t_min, time_it_t const& t_max,
                                      dcomplex c) {
 double dt = *t_max - *t_min;
 if(!inst_h.is_static) diag(*t_min, dt);
 auto & psi = st.amplitudes();
 psi = conj(eig.second) * psi;
 for(long i : range(psi.size())) psi(i) *= std::exp(c * dt * eig.first(i));
 psi = eig.second.transpose() * psi;
}

//////////////////////
/// prop_lanczos_t ///
//////////////////////

inline prop_lanczos_t::prop_lanczos_t(inst_h_t && inst_h_, sub_hilbert_space const& sp,
                                      double gs_energy_convergence, int max_krylov_dim) :
 inst_h(std::move(inst_h_)),
 lw(inst_h, gs_energy_convergence > 0 ? gs_energy_convergence : 1e-8,
            std::min((max_krylov_dim > 0 ? max_krylov_dim : 20), sp.size())),
 lanczos_exp(lw.max_krylov_dim, lw.max_krylov_dim) {
}

inline void prop_lanczos_t::operator()(state_on_subspace_t & st,
                                       time_it_t const& t_min, time_it_t const& t_max,
                                       dcomplex c) {
 inst_h.t = *t_min;
 inst_h.dt = *t_max - *t_min;

 // Construct the Krylov basis
 double norm = std::sqrt(lw.checked_real(dot_product(st, st)));
 st /= norm;
 lw(st);

 // Construct the propagation exponent
 vector_view<double> eigenvalues = lw.values();
 auto all = range(eigenvalues.size());

 for(int n : all)
  lanczos_exp(n, all) = std::exp(c * inst_h.dt * eigenvalues(n)) * lw.vectors()(n, all);
 lanczos_exp(all, all) = lw.vectors().transpose() * lanczos_exp(all, all);

 // Propagate
 lw.krylov_2_fock(norm * lanczos_exp(all, 0), st);
}

//////////////////
/// propagator ///
//////////////////

propagator::propagator(op_on_subspace_t const& h, sub_hilbert_space const& sp,
                       double hbar, bool is_static_h, h_interpolation h_interpol,
                       long lanczos_min_matrix_size,
                       double gs_energy_convergence,
                       int max_krylov_dim
                      ) :
 is_static_h(is_static_h), h_coeff(-1_j / hbar),
 ed_solver(sp.size() == 1 ? OneD : (sp.size() < lanczos_min_matrix_size ? LAPACK : Lanczos)) {
 switch(ed_solver) {
  case OneD:
   prop_impl = std::make_shared<prop_1d_t>(inst_h_t(h, is_static_h, h_interpol), sp);
   break;
  case LAPACK:
   prop_impl = std::make_shared<prop_lapack_t>(inst_h_t(h, is_static_h, h_interpol), sp);
   break;
  case Lanczos:
   prop_impl = std::make_shared<prop_lanczos_t>(inst_h_t(h, is_static_h, h_interpol), sp,
                                                gs_energy_convergence, max_krylov_dim);
   break;
 }
}

inline void propagator::propagate(state_on_subspace_t & st,
                                  time_it_t const& t_min,
                                  time_it_t const& t_max,
                                  dcomplex c) const {
 switch(ed_solver) {
  case OneD:    (*static_cast<prop_1d_t*>(prop_impl.get()))(st, t_min, t_max, c); break;
  case LAPACK:  (*static_cast<prop_lapack_t*>(prop_impl.get()))(st, t_min, t_max, c); break;
  case Lanczos: (*static_cast<prop_lanczos_t*>(prop_impl.get()))(st, t_min, t_max, c); break;
 }
}

void propagator::operator()(state_on_subspace_t & st,
                            time_it_t const& t_start,
                            time_it_t const& t_end) const {
 if(t_end == t_start) return; // No propagation needed

 bool forward = *t_end > *t_start;
 dcomplex c = (forward ? 1 : -1) * h_coeff;

 if(is_static_h && ed_solver != Lanczos) {
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

}
