/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2026, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <limits>
#include <utility>

#include <h5/h5.hpp>
#include <mpi/mpi.hpp>

#include "types.hpp"

namespace realevol {

class init_state;

/// Apply a generating operator to the vacuum state to make a pure initial state
init_state make_pure_init_state(static_operator_t const& generator,
                                fundamental_operator_set const& fops,
                                std::map<operators::indices_t, int> const& bits_per_boson = {});

struct eq_solver_parameters_t {

 /// Verbosity level for the equilibrium solver
 int verbosity = 0;

 /// Discard states with relative statistical weight below this threshold
 double min_rel_weight = std::numeric_limits<double>::epsilon();

 /// Call ARPACK to diagonalize matrices of this size or bigger (cannot be < 4)
 int arpack_min_matrix_size = 101;

 /// Eigenvalue convergence tolerance for ARPACK
 double arpack_tolerance = 0;

 /// ARPACK parameter NCV (number of Lanczos vectors) for each invariant subspace
 std::map<long,int> arpack_ncv = std::map<long,int>({});
};

/// Make an equilibrium initial state at a given temperature (zero temperature is valid)
init_state make_equilibrium_init_state(static_operator_t const& h,
                                       fundamental_operator_set const& fops,
                                       double temperature,
                                       eq_solver_parameters_t const& params,
                                       std::map<operators::indices_t, int> const& bits_per_boson = {},
                                       mpi::communicator comm = {});

/// Initial state, including information about the Hilbert space structure
class init_state {

 /// Fundamental operator set used to construct this state
 fundamental_operator_set fops;
 /// Full Hilbert space
 class hilbert_space full_hs;
 /// Invariant subspaces with relevant contributions to this state
 std::vector<sub_hilbert_space> sub_hilbert_spaces;

 // Constructor to be used only by the factory functions
 init_state(fundamental_operator_set const& fops, std::map<operators::indices_t, int> const& bits_per_boson = {}) :
  fops(fops) {
  if(fops.size(Boson) != bits_per_boson.size()) TRIQS_RUNTIME_ERROR
   << "bits_per_boson has a wrong number of elements, must be " << fops.size(Boson);
  std::vector<int> v(bits_per_boson.size());
  for(auto const& b : bits_per_boson) {
   if(!fops.has_indices(b.first, Boson)) TRIQS_RUNTIME_ERROR
    << "Wrong indices " << b.first << " in bits_per_boson";
   v[fops.pos(b.first, Boson)] = b.second;
  }
  full_hs = {fops, v};
 }

public:

 /// Constructor
 init_state() = default;

 // Do not allow copying
 init_state(init_state const&) = delete;
 init_state & operator=(init_state const&) = delete;

 /// Move-constructor
 init_state(init_state &&) noexcept = default;

 /// Move-assignment operator
 init_state & operator=(init_state && ist) noexcept {
  using std::swap;
  swap(fops,               ist.fops);
  swap(full_hs,            ist.full_hs);
  swap(sub_hilbert_spaces, ist.sub_hilbert_spaces);
  swap(weighted_states,    ist.weighted_states);
  return *this;
 }

 struct weighted_state_t {
  /// State, in one of the subspaces
  state_on_subspace_t state;
  /// Statistical weight
  double weight;

  /// Constructor
  weighted_state_t(state_on_subspace_t && state, double weight) :
   state(std::move(state)), weight(weight) {}
 };

 /// Constant iterator over weighted pure states
 using const_iterator = std::vector<weighted_state_t>::const_iterator;

 [[nodiscard]] inline const_iterator begin() const noexcept { return weighted_states.cbegin(); }
 [[nodiscard]] inline const_iterator cbegin() const noexcept { return weighted_states.cbegin(); }
 [[nodiscard]] inline const_iterator end() const noexcept { return weighted_states.cend(); }
 [[nodiscard]] inline const_iterator cend() const noexcept { return weighted_states.cend(); }

 /// Get fundamental operator set
 [[nodiscard]] fundamental_operator_set const& get_fops() const { return fops; }

 /// Get full Hilbert space
 [[nodiscard]] class hilbert_space const& get_full_hs() const { return full_hs; }

 /// Get invariant subspaces
 [[nodiscard]] std::vector<sub_hilbert_space> const& get_sub_hilbert_spaces() const {
  return sub_hilbert_spaces;
 }

 /// Number of weighted pure states in this initial state
 [[nodiscard]] int size() const { return weighted_states.size(); }

 /// Get all weighted states
 [[nodiscard]] std::vector<weighted_state_t> const& get_weighted_states() const { return weighted_states; }

 /// Stream output
 friend std::ostream & operator<<(std::ostream & os, init_state const& st);

 /// HDF5 read/write
 static std::string hdf5_format() { return "InitState"; }
 friend void h5_write(h5::group gr, std::string const &name, init_state const& st);
 friend void h5_read(h5::group gr, std::string const &name, init_state & st);

private:

 /// List of pure states
 std::vector<weighted_state_t> weighted_states;

 friend init_state make_pure_init_state(static_operator_t const&,
                                        fundamental_operator_set const&,
                                        std::map<operators::indices_t, int> const&);

 friend init_state make_equilibrium_init_state(static_operator_t const&,
                                               fundamental_operator_set const&,
                                               double,
                                               eq_solver_parameters_t const&,
                                               std::map<operators::indices_t, int> const&,
                                               mpi::communicator);

};

}
