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
#pragma once

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <utility>
#include <triqs/h5.hpp>

#include "common.hpp"

namespace realevol {

class init_state;

/// Apply a generating operator to the vacuum state to make a pure initial state
init_state make_pure_init_state(operator_t const& generator,
                                fundamental_operator_set const& fops,
                                std::map<operators::indices_t, int> const& bits_per_boson = {});

/// FIXME
init_state make_zerotemp_init_state();

/// FIXME
init_state make_thermal_init_state();

// std::vector<init_state> init_statehermal(operator_t const& h0, hilbert_space_structure & hss,
//                                              double beta,
//                                              double temperature_cutoff = std::numeric_limits<double>::epsilon());

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
  TRIQS_ASSERT(fops.size(Boson) == bits_per_boson.size());
  std::vector<int> v(bits_per_boson.size());
  for(auto const& b : bits_per_boson) v[fops.pos(b.first,Boson)] = b.second;
  full_hs = {fops, v};
 }

public:

 /// Constructor
 init_state() = default;

 struct weighted_state_t {
  /// State, in one of the subspaces
  state_on_subspace_t state;
  /// Statistical weight
  double weight;

  weighted_state_t() = default;

  /// Constructor
  weighted_state_t(state_on_subspace_t && state, double weight) :
   state(std::move(state)), weight(weight) {}
 };

 /// Constant iterator over weighted pure states
 using const_iterator = std::vector<weighted_state_t>::const_iterator;

 inline const_iterator begin() const noexcept { return weighted_states.cbegin(); }
 inline const_iterator cbegin() const noexcept { return weighted_states.cbegin(); }
 inline const_iterator end() const noexcept { return weighted_states.cend(); }
 inline const_iterator cend() const noexcept { return weighted_states.cend(); }

 /// Get fundamental operator set
 fundamental_operator_set const& get_fops() const {return fops; }

 /// Get full Hilbert space
 class hilbert_space const& get_full_hs() const { return full_hs; }

 /// Get invariant subspaces
 std::vector<sub_hilbert_space> const& get_sub_hilbert_spaces() const {
  return sub_hilbert_spaces;
 }

 /// Number of weighted pure states in this initial state
 int size() const { return weighted_states.size(); }

 /// Get all weighted states
 std::vector<weighted_state_t> const& get_weighted_states() const { return weighted_states; }

 /// Stream output
 friend std::ostream & operator<<(std::ostream & os, init_state const& st);

 /// HDF5 read/write
 friend std::string get_triqs_hdf5_data_scheme(init_state const&) { return "InitState"; }
 friend void h5_write(h5::group gr, std::string const &name, init_state const& st);
 friend void h5_read(h5::group gr, std::string const &name, init_state & st);

private:

 /// List of pure states
 std::vector<weighted_state_t> weighted_states;

 friend init_state make_pure_init_state(operator_t const&,
                                        fundamental_operator_set const&,
                                        std::map<operators::indices_t, int> const&);

};

}
