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

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

#include <cmath>

#include <triqs/arrays/linalg/eigenelements.hpp>

#include "time_expr.hpp"
#include "time_interp.hpp"
#include "propagator.hpp"

namespace realevol {

using namespace triqs::arrays;

////////////////
/// inst_h_t ///
////////////////

template<typename HScalarType>
inline inst_h_t<HScalarType>::inst_h_t(op_on_subspace_t<HScalarType> const& h,
                                       bool is_static,
                                       h_interpolation h_interpol) :
 is_static(is_static), h(h), h_interpol(is_static ? Rectangle : h_interpol), t(0), dt(0) {
}

template<typename HScalarType>
inline state_on_subspace_t inst_h_t<HScalarType>::operator()(state_on_subspace_t const& st) const {
 switch(h_interpol) {
  case Rectangle: return h(st, t+dt/2);
  case Trapezoid: return (h(st, t) + h(st, t+dt)) / 2.0;
  case Simpson:   return (h(st, t) + 4*h(st, t+dt/2) + h(st, t+dt)) / 6.0;
 }
}

/////////////////
/// prop_1d_t ///
/////////////////

template<typename HScalarType>
inline prop_1d_t<HScalarType>::prop_1d_t(inst_h_t<HScalarType> && inst_h_, sub_hilbert_space const& sp) :
 inst_h(std::move(inst_h_)), tmp_st(sp) {
 tmp_st(0) = 1;
 if(inst_h.is_static) diag(0,0);
}

template<typename HScalarType>
inline void prop_1d_t<HScalarType>::diag(double t, double dt) {
 inst_h.t = t;
 inst_h.dt = dt;
 eigenvalue = inst_h(tmp_st)(0).real();
}

template<typename HScalarType>
inline void prop_1d_t<HScalarType>::operator()(state_on_subspace_t & st, double t_min, double t_max, dcomplex c) {
 double dt = t_max - t_min;
 if(!inst_h.is_static) diag(t_min, dt);
 st *= std::exp(c * dt * eigenvalue);
}

/////////////////////
/// prop_lapack_t ///
/////////////////////

template<typename HScalarType>
inline prop_lapack_t<HScalarType>::prop_lapack_t(inst_h_t<HScalarType> && inst_h_, sub_hilbert_space const& sp) :
 inst_h(std::move(inst_h_)), from_st(sp), to_st(sp), workspace(sp.size(), sp.size()) {
 if(inst_h.is_static) diag(0,0);
}

template<typename HScalarType>
inline void prop_lapack_t<HScalarType>::diag(double t, double dt) {
 inst_h.t = t;
 inst_h.dt = dt;
 for(int i : range(from_st.size())) {
  from_st(i) = 1.0;
  to_st = inst_h(from_st);
  workspace(range(),i) = to_st.amplitudes();
  from_st(i) = 0;
 }
 eig = linalg::eigenelements_in_place(&workspace);
}

template<typename HScalarType>
inline void prop_lapack_t<HScalarType>::operator()(state_on_subspace_t & st, double t_min, double t_max, dcomplex c) {
 double dt = t_max - t_min;
 if(!inst_h.is_static) diag(t_min, dt);
 auto & psi = st.amplitudes();
 psi = conj(eig.second) * psi;
 for(int i : range(psi.size())) psi(i) *= std::exp(c * dt * eig.first(i));
 psi = eig.second.transpose() * psi;
}

//////////////////////
/// prop_lanczos_t ///
//////////////////////

template<typename HScalarType>
inline prop_lanczos_t<HScalarType>::prop_lanczos_t(inst_h_t<HScalarType> && inst_h_, sub_hilbert_space const& sp,
                                                   double gs_energy_convergence, int max_krylov_dim) :
 inst_h(std::move(inst_h_)),
 lw(inst_h, gs_energy_convergence > 0 ? gs_energy_convergence : 1e-8,
            std::min((max_krylov_dim > 0 ? max_krylov_dim : 20), sp.size())),
 krylov_coeffs(lw.max_krylov_dim) {
}

template<typename HScalarType>
inline void prop_lanczos_t<HScalarType>::operator()(state_on_subspace_t & st, double t_min, double t_max, dcomplex c) {
 inst_h.t = t_min;
 inst_h.dt = t_max - t_min;

 // Construct the Krylov basis
 double norm = std::sqrt(lw.checked_real(dot_product(st, st)));
 st /= norm;
 lw(st);

 // Construct the propagation exponent
 vector_view<double> eigenvalues = lw.values();
 auto all = range(eigenvalues.size());

 for(int n : all)
  krylov_coeffs(n) = std::exp(c * inst_h.dt * eigenvalues(n)) * lw.vectors()(n, 0);
 krylov_coeffs(all) = lw.vectors().transpose() * krylov_coeffs(all);

 // Propagate
 lw.krylov_2_fock(norm * krylov_coeffs(all), st);
}

//////////////////
/// propagator ///
//////////////////

template<typename HScalarType>
propagator<HScalarType>::propagator(op_on_subspace_t<HScalarType> const& h, sub_hilbert_space const& sp,
                                    gf_mesh<retime> const& t_mesh,
                                    double hbar, bool is_static_h, h_interpolation h_interpol,
                                    long lanczos_min_matrix_size,
                                    double gs_energy_convergence,
                                    int max_krylov_dim
                                    ) :
 is_static_h(is_static_h), h_coeff(-1i / hbar), t_mesh(t_mesh),
 ed_solver(sp.size() == 1 ? OneD : (sp.size() < lanczos_min_matrix_size ? LAPACK : Lanczos)) {
 switch(ed_solver) {
  case OneD:
   prop_impl = std::make_shared<prop_1d_t<HScalarType>>(inst_h_t<HScalarType>(h, is_static_h, h_interpol), sp);
   break;
  case LAPACK:
   prop_impl = std::make_shared<prop_lapack_t<HScalarType>>(inst_h_t<HScalarType>(h, is_static_h, h_interpol), sp);
   break;
  case Lanczos:
   prop_impl = std::make_shared<prop_lanczos_t<HScalarType>>(inst_h_t<HScalarType>(h, is_static_h, h_interpol), sp,
                                                             gs_energy_convergence, max_krylov_dim);
   break;
 }
}


template<typename HScalarType>
inline void propagator<HScalarType>::propagate(state_on_subspace_t & st, int t_min_index, int t_max_index, dcomplex c) const {
 switch(ed_solver) {
  case OneD:    (*static_cast<prop_1d_t<HScalarType>*>(prop_impl.get()))(st, t_mesh[t_min_index], t_mesh[t_max_index], c); break;
  case LAPACK:  (*static_cast<prop_lapack_t<HScalarType>*>(prop_impl.get()))(st, t_mesh[t_min_index], t_mesh[t_max_index], c); break;
  case Lanczos: (*static_cast<prop_lanczos_t<HScalarType>*>(prop_impl.get()))(st, t_mesh[t_min_index], t_mesh[t_max_index], c); break;
 }
}

template<typename HScalarType>
void propagator<HScalarType>::operator()(state_on_subspace_t & st, int t_start_index, int t_end_index) const {
 if(t_end_index == t_start_index) return; // No propagation needed

 if(is_static_h && ed_solver != Lanczos) {
  propagate(st, t_start_index, t_end_index, h_coeff); // Forward or backward propagation (static)
 } else {
  bool forward = t_end_index > t_start_index;
  if(forward) { // Forward propagation
   int t_next_index = t_start_index; ++t_next_index;
   for(int t_index = t_start_index; t_index != t_end_index; ++t_index, ++t_next_index)
    propagate(st, t_index, t_next_index, h_coeff);
  } else { // Backward propagation
   int t_next_index = t_start_index; --t_next_index;
   for(auto t_index = t_start_index; t_index != t_end_index; --t_index, --t_next_index)
    propagate(st, t_index, t_next_index, h_coeff);
  }
 }
}

template class propagator<time_expr>;
template class propagator<time_interp>;

}
