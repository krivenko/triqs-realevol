/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013-2020, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include <algorithm>
#include <cmath>

#include <boost/operators.hpp>
#include <triqs/utility/numeric_ops.hpp>
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>

#include "hilbert_space.hpp"
#include "state.hpp"

namespace realevol {
namespace hilbert_space {

/// Many-body state as a view of a `triqs::arrays::vector` object.
/**
  @tparam HilbertSpace Hilbert space type, one of [[hilbert_space]] and [[sub_hilbert_space]]
  @tparam ScalarType Amplitude type, normally `double` or `std::complex<double>`
  @include triqs/hilbert_space/state_view.hpp
 */
template <typename HilbertSpace, typename ScalarType>
class state_view : boost::additive<state_view<HilbertSpace, ScalarType>>,
                   boost::multiplicative<state_view<HilbertSpace, ScalarType>, ScalarType> {

 const HilbertSpace* hs_p;
 using amplitude_t = triqs::arrays::vector_view<ScalarType>;
 amplitude_t ampli;

 public:
 /// Accessor to `ScalarType` template parameter
 using value_type = ScalarType;
 /// Accessor to `HilbertSpace` template parameter
 using hilbert_space_t = HilbertSpace;

 /// Construct a new state_view object
 /**
   @param hs Hilbert space the new state_view belongs to
  */
 state_view(amplitude_t ampli, HilbertSpace const& hs) : hs_p(&hs), ampli(ampli) {
  TRIQS_ASSERT(hs.size() == ampli.size());
 }

 /// Copy-constructor
 state_view(state_view const&) = default;

 /// Return the dimension of the associated Hilbert space
 /**
   @return Dimension of the associated Hilbert space
  */
 int size() const { return hs_p->size(); }

 /// Access to individual amplitudes
 /**
   @param i index of the requested amplitude
   @return Reference to the requested amplitude
  */
 value_type& operator()(int i) { return ampli[i]; }
 /// Access to individual amplitudes
 /**
   @param i index of the requested amplitude
   @return Constant reference to the requested amplitude
  */
 value_type const& operator()(int i) const { return ampli[i]; }

 /// In-place addition of another state_view
 /**
   @param s2 Another [[state_view]] object to add
   @return Reference to this state_view
  */
 state_view& operator+=(state_view const& s2) {
  ampli += s2.ampli;
  return *this;
 }

 /// In-place subtraction of another state_view
 /**
   @param s2 Another [[state_view]] object to add
   @return Reference to this state_view
  */
 state_view& operator-=(state_view const& s2) {
  ampli -= s2.ampli;
  return *this;
 }

 /// In-place multiplication by a scalar
 /**
   @param x Multiplier
   @return Reference to this state_view
  */
 state_view& operator*=(value_type x) {
  ampli *= x;
  return *this;
 }

 /// In-place division by a scalar
 /**
   @param x Divisor
   @return Reference to this state_view
  */
 state_view& operator/=(value_type x) {
  ampli /= x;
  return *this;
 }

 /// Set all amplitdes to zero
 void zero() { ampli() = value_type{}; }

 /// Calculate scalar product of two states
 /**
   @param sv1 First state_view to multiply
   @param sv2 Second state_view to multiply
   @return Value of the scalar product
  */
 friend value_type dot_product(state_view const& sv1, state_view const& sv2) {
  return dotc(sv1.ampli, sv2.ampli);
 }

 /// Apply a callable object to all amplitudes of a state_view
 /**
  The callable must take two arguments, 1) index of the basis Fock state in the associated Hilbert space, and 2) the corresponding amplitude.

  @tparam Lambda Type of the callable object
  @param sv [[state_view]] object
  @param l Callable object
  */
 template<typename Lambda>
 friend void foreach(state_view const& sv, Lambda l) {
  const auto L = sv.size();
  for (size_t i = 0; i < L; ++i) l(i, sv(i));
 }

 //
 // Additions to StateVector concept
 //

 /// Direct access to the storage container (`triqs::arrays::vector_view`)
 /**
   @return Constant reference to the storage container
  */
 amplitude_t const& amplitudes() const { return ampli; }
 /// Direct access to the storage container (`triqs::arrays::vector`)
 /**
   @return Reference to the storage container
  */
 amplitude_t& amplitudes() { return ampli; }

 /// Return a constant reference to the associated Hilbert space
 /**
   @return Constant reference to the Hilbert space
  */
 HilbertSpace const& get_hilbert() const { return *hs_p; }
 /// Reset the associated Hilbert space
 /**
   @param new_hs Constant reference to the new Hilbert space
  */
 void set_hilbert(HilbertSpace const& new_hs) { hs_p = &new_hs; }

};

// Print state_view
template <typename HilbertSpace, typename ScalarType>
std::ostream& operator<<(std::ostream& os, state_view<HilbertSpace, ScalarType> const& sv) {
 bool something_written = false;
 auto const& hs = sv.get_hilbert();

 using value_type = typename state_view<HilbertSpace, ScalarType>::value_type;
 foreach(sv, [&os,hs,&something_written](int i, value_type ampl){
  using triqs::utility::is_zero;
  if (!is_zero(ampl)){
   os << " +(" << ampl << ")" << "|" << hs.get_fock_state(i) << ">";
   something_written = true;
  }
 });

 if (!something_written) os << 0;
 return os;
}

/// Make a state object equivalent to a given state_view with all amplitudes set to 0
/**
  @param st State view object
  @tparam HilbertSpace Hilbert space type, one of [[hilbert_space]] and [[sub_hilbert_space]]
  @tparam ScalarType Amplitude type, normally `double` or `std::complex<double>`
  @include triqs/hilbert_space/view_state.hpp
 */
template <typename HilbertSpace, typename ScalarType, bool BasedOnMap>
state<HilbertSpace, ScalarType, false> make_zero_state(state_view<HilbertSpace, ScalarType> const& st) {
 return {st.get_hilbert()};
}

}}
