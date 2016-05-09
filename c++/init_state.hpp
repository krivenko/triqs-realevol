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

#include <vector>
#include <utility>
#include <algorithm>
#include <triqs/hilbert_space/state.hpp>

#include "hs_structure.hpp"

namespace realevol {

using state_on_subspace_t = state<sub_hilbert_space,dcomplex,false>;
using init_states_t = std::vector<std::pair<state_on_subspace_t,double>>;

bool check_operator_static(operator_t const& op) {
 for(auto const& m : op)
  if(!is_constant(m.coef)) return false;
 return true;
}

init_states_t init_states_pure(operator_t const& generator, hilbert_space_structure & hss) {

 // Check if the generating operator is static
 if(!check_operator_static(generator))
  TRIQS_RUNTIME_ERROR << "Generating operator must be time-independent!";

 init_states_t states;
 std::vector<state_on_subspace_t*> state_part_ptrs(hss.sub_hilbert_spaces.size(), nullptr);

 state_on_space_t vacuum(hss.full_hs);
 vacuum(fock_state_t(0)) = 1;

 op_on_space_t imp_gen(generator, hss.fops, hss.full_hs);
 auto init_state = imp_gen(vacuum, 0);

 foreach(init_state, [&states,&state_part_ptrs,&hss](fock_state_t f, dcomplex a){
  auto spn = hss.partition.lookup_basis_state(f);
  auto & ptr = state_part_ptrs[spn];
  auto const& sp = hss.sub_hilbert_spaces[spn];
  if(ptr == nullptr) {
   states.emplace_back(sp, 1.0);
   ptr = &states.back().first;
  }
  (*ptr)(sp.get_state_index(f)) = a;
 });

 return states;
}

}
