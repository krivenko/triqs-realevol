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
#include <triqs/utility/first_include.hpp>

#include "init_state.hpp"

namespace realevol {

bool check_operator_static(operator_t const& op) {
 for(auto const& m : op)
  if(!is_constant(m.coef)) return false;
 return true;
}

std::vector<init_state_t> init_state_pure(operator_t const& generator, hilbert_space_structure & hss) {

 // Check if the generating operator is static
 if(!check_operator_static(generator))
  TRIQS_RUNTIME_ERROR << "Generating operator must be time-independent!";

 std::vector<init_state_t> states = {{{}, 1.0}};
 auto & parts = states[0].parts;

 state_on_space_t vacuum(hss.full_hs);
 vacuum(fock_state_t(0)) = 1;

 op_on_space_t imp_gen(generator, hss.fops, hss.full_hs);
 auto init_state = imp_gen(vacuum, 0);

 foreach(init_state, [&parts,&hss](fock_state_t f, dcomplex a){
  auto spn = hss.partition.lookup_basis_state(f);
  auto const& sp = hss.sub_hilbert_spaces[spn];
  auto ins = parts.emplace(spn, sp);
  (ins.first->second)(sp.get_state_index(f)) = a;
 });

 return states;
}

}
