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

inline prop_1d_t::prop_1d_t(inst_h_t && inst_h, sub_hilbert_space const& sp) :
 inst_h(inst_h), tmp_st(sp) {
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

inline prop_lapack_t::prop_lapack_t(inst_h_t && inst_h, sub_hilbert_space const& sp) :
 inst_h(inst_h), from_st(sp), to_st(sp), workspace(sp.size(), sp.size()) {
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
/// diag_lanczos_t ///
//////////////////////
inline prop_lanczos_t::prop_lanczos_t(inst_h_t && inst_h,
                                      double gs_energy_convergence, int max_krylov_dim) :
 inst_h(inst_h), lw(inst_h, gs_energy_convergence, max_krylov_dim) {
}

inline void prop_lanczos_t::operator()(state_on_subspace_t & st,
                                       time_it_t const& t_min, time_it_t const& t_max,
                                       dcomplex c) {

// TODO
/*
-        auto next = first; ++next;
-
-        // Value of the independent variable
-        var_t x;
-
-        // RHS as a function of the solution only
-        using std::placeholders::_1;
-        auto rhs_of_sol_only = std::bind(hermitian_rhs, _1, std::cref(x));
-
-        // Lanczos worker object
-        lanczos_worker<decltype(rhs_of_sol_only),value_t> lw(rhs_of_sol_only,gs_energy_convergence);
-
-        // Dimension of the problem
-        std::size_t N = (first->value).size();
-
-        triqs::arrays::matrix<scalar_t> lanczos_exp(N,N);
-
-        for(;next != last; ++first, ++next){
-            // Value of the independent variable
-            x = 0.5*(first->mesh_point + next->mesh_point);
-
-            // Solution at current point
-            value_t const& U(first->value);
-
-            // Construct the Krylov basis
-            auto norm = std::sqrt(dot_product(U,U));
-            lw(U/norm);
-
-            // Step of the mesh
-            var_t step = next->mesh_point - first->mesh_point;
-
-            // Construct the propagation exponent
-            auto eigenvalues = lw.values();
-            std::size_t krylov_dim = eigenvalues.size();
-            auto all = range(0,krylov_dim);
-
-            for (std::size_t n = 0; n < krylov_dim; ++n)
-                lanczos_exp(n,all) = exp(step * rhs_prefactor * eigenvalues(n)) * lw.vectors()(n,all);
-            lanczos_exp(all,all) = lw.vectors().transpose() * lanczos_exp(all,all);
-
-            // Propagate
-            auto krylov_coeffs = norm * lanczos_exp(all, 0);
-            next->value = lw.krylov_2_fock(krylov_coeffs);
-        }
-    }
*/
}

//////////////////
/// propagator ///
//////////////////

std::shared_ptr<void> propagator::make_prop_impl(op_on_subspace_t const& h, sub_hilbert_space const& sp,
                                                 bool is_static_h, h_interpolation h_interpol) {
 inst_h_t inst_h {h, h_interpol, .0, .0, is_static_h};
 switch(ed_solver) {
  case OneD: return std::make_shared<prop_1d_t>(std::move(inst_h), sp);
  case LAPACK: return std::make_shared<prop_lapack_t>(std::move(inst_h), sp);
  case Lanczos: return std::make_unique<prop_lanczos_t>(std::move(inst_h), 1e-8, 20 /* FIXME */);
 }
}

propagator::propagator(op_on_subspace_t const& h, sub_hilbert_space const& sp,
                       double hbar, bool is_static_h, h_interpolation h_interpol,
                       long lanczos_min_matrix_size) :
 is_static_h(is_static_h), h_coeff(-1_j / hbar),
 ed_solver(sp.size() == 1 ? OneD : (sp.size() < lanczos_min_matrix_size ? LAPACK : Lanczos)),
 prop_impl(make_prop_impl(h, sp, is_static_h, h_interpol)) {
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

inline void propagator::propagate(state_on_subspace_t & st,
                                  time_it_t const& t_min,
                                  time_it_t const& t_max,
                                  dcomplex c) const {
 switch(ed_solver) {
  case OneD:    (*static_cast<prop_1d_t*>(prop_impl.get()))(st, t_min, t_max, c);
  case LAPACK:  (*static_cast<prop_lapack_t*>(prop_impl.get()))(st, t_min, t_max, c);
  case Lanczos: (*static_cast<prop_lanczos_t*>(prop_impl.get()))(st, t_min, t_max, c);
 }
}

}
