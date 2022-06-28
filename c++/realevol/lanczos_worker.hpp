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
#pragma once

#include <complex>
#include <type_traits>
#include <vector>

#include <nda/nda.hpp>

#include <triqs/utility/is_complex.hpp>
#include <triqs/utility/numeric_ops.hpp>
#include <triqs/utility/exceptions.hpp>

#include <triqs/tridiag_worker.hpp>

using namespace nda;
using nda::lapack::tridiag_worker;
using triqs::is_complex;
using triqs::utility::is_zero;

namespace realevol {

template <typename OperatorType, typename StateType> struct lanczos_worker {

 using scalar_t = typename StateType::value_type;
 using real_scalar_t = decltype(std::real(scalar_t{}));

 OperatorType const& h;

 // Krylov basis states
 // h \approx V * T * V^+
 // Elements of 'basisstates' are columns of V
 std::vector<StateType> basisstates;

 // The tridiagonal matrix T has the following form
 // | alpha[0]    beta[0]     0       ...     |
 // | beta[0]     alpha[1]    beta[1]     ... |
 // | 0            beta[1]    alpha[2]    ... |
 // |             ...                         |
 std::vector<real_scalar_t> alpha; // diagonal matrix elements
 std::vector<real_scalar_t> beta;  // superdiagonal matrix elements

 // Temporaries
 StateType res_vector;

 // Convergence threshold for the GS energy
 real_scalar_t gs_energy_convergence;

 // Maximal dimension of the Krylov space
 int max_krylov_dim;

 // Tridiagonal matrix diagonalizer
 tridiag_worker<false> tdw;

 // For complex numbers: extract the real part and make sure that the imaginary part is negligible
 // For real numbers: returns the argument
 template<bool IsComplex = is_complex<scalar_t>::value>
 inline real_scalar_t checked_real(scalar_t const& x,
                                   std::enable_if_t<IsComplex,void*> = 0) {
  TRIQS_ASSERT(is_zero(std::imag(x)));
  return std::real(x);
 }
 template<bool IsComplex = is_complex<scalar_t>::value>
 inline real_scalar_t checked_real(scalar_t x,
                                   std::enable_if_t<!IsComplex,void*> = 0) {
  return x;
 }

 // Returns the only matrix element of the 1x1 Krylov-projected matrix
 real_scalar_t first_iteration(StateType const& initial_state) {
  basisstates.push_back(initial_state);
  res_vector = h(basisstates.back());
  alpha.push_back(checked_real(dot_product(initial_state, res_vector)));
  res_vector -= alpha.back() * initial_state;
  return alpha.back();
 }

 // Calculates the next state in Krylov's basis.
 // Returns false if the previous state was an eigenstate of h
 bool advance() {
  real_scalar_t new_beta = std::sqrt(checked_real(dot_product(res_vector, res_vector)));
  // We don't really want to divide by zero
  if(is_zero(new_beta, gs_energy_convergence)) return false;
  beta.push_back(new_beta);
  basisstates.push_back(res_vector / new_beta);
  res_vector = h(basisstates.back());
  alpha.push_back(checked_real(dot_product(basisstates.back(), res_vector)));
  res_vector -= alpha.back() * basisstates.back();
  res_vector -= beta.back() * basisstates[basisstates.size() - 2];
  return true;
 }

public:

 using state_type = StateType;

 lanczos_worker(OperatorType const& h, real_scalar_t gs_energy_convergence, int max_krylov_dim = 20) :
  h(h), gs_energy_convergence(gs_energy_convergence), max_krylov_dim(max_krylov_dim),
  tdw(max_krylov_dim) {
  TRIQS_ASSERT(max_krylov_dim >= 1);
  alpha.reserve(max_krylov_dim);
  beta.reserve(max_krylov_dim - 1);
  basisstates.reserve(max_krylov_dim);
 }
 lanczos_worker(lanczos_worker const&) = default;
 lanczos_worker& operator=(lanczos_worker const&) = delete;

 // initial_state MUST be of norm 1
 void operator()(StateType const& initial_state) {
  reset();

  // First iteration
  real_scalar_t gs_energy = first_iteration(initial_state);

  tdw(alpha, beta); // FIXME: no need to call this... but otherwise tdw.values() can not be called

  while (advance()) {
   tdw(alpha, beta);
   if(is_zero(tdw.values()[0] - gs_energy, gs_energy_convergence)) break;
   if(tdw.values().size() == max_krylov_dim) {
    //std::cerr << "lanczos_worker: maximal dimension of the Krylov space ("
    //          << max_krylov_dim << ") has been reached"  << std::endl;
    break;
   }
   gs_energy = tdw.values()[0];
  }
 }

 // Access eigenvalues and eigenvectors of the Krylov-projected operator
 vector_view<double> values() const { return tdw.values(); }
 matrix_view<double> vectors() const { return tdw.vectors(); }

 void reset() {
  alpha.clear();
  beta.clear();
  basisstates.clear();
 }

 template <typename KrylovCoeffs> void krylov_2_fock(KrylovCoeffs const& phi, StateType & st) {
  st.zero();
  for (int i = 0; i < phi.size(); ++i) st += phi(i) * basisstates[i];
 }

 template <typename KrylovCoeffs> StateType krylov_2_fock(KrylovCoeffs const& phi) {
  state_type st = make_zero_state(res_vector);
  for (int i = 0; i < phi.size(); ++i) st += phi(i) * basisstates[i];
  return st;
 }
};

}

