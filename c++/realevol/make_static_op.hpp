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

#include "types.hpp"

namespace realevol {

template<typename OperatorType>
static_operator_t make_static_op(OperatorType const& op) {
 static_operator_t static_op;
 for(auto const& m : op) {
  if(!is_constant(m.coef)) TRIQS_RUNTIME_ERROR << "Operator must be time-independent";
  auto new_coef = m.coef(0);
  static_operator_t new_monomial(new_coef);
  for(auto const& c : m.monomial)
   new_monomial *= static_operator_t::make_canonical(c.stat, c.dagger, c.indices);
  static_op += new_monomial;
 }
 return static_op;
}

} // namespace realevol
