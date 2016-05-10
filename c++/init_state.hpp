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

#include <map>
#include <limits>

#include "hs_structure.hpp"

namespace realevol {

// Initial state
struct init_state_t {
 // subspace index -> projection of this state onto the subspace
 std::map<int, state_on_subspace_t> parts;
 // Statistical weight of this initial state
 double weight;
};

std::vector<init_state_t> init_state_pure(operator_t const& generator, hilbert_space_structure & hss);
std::vector<init_state_t> init_state_thermal(operator_t const& h0, hilbert_space_structure & hss,
                                             double beta,
                                             double temperature_cutoff = std::numeric_limits<double>::epsilon());

}
