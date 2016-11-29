/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2015 by O. Parcollet
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
#include "./many_body_operator.hpp"
#include <triqs/h5.hpp>
#include <triqs/h5/base.hpp>

namespace realevol {
namespace operators {

std::ostream& operator<<(std::ostream& os, canonical_ops_t const& op) {
 os << (op.stat == hilbert_space::statistic_enum::Fermion ? "C" : "A");
 if (op.dagger) os << "^+";
 os << "(";
 int u = 0;
 for (auto const& i : op.indices) {
  if (u++) os << ",";
  os << i;
 }
 return os << ")";
}

bool operator<(monomial_t const& m1, monomial_t const& m2) {
 return m1.size() != m2.size() ? m1.size() < m2.size()
                               : std::lexicographical_compare(m1.begin(), m1.end(), m2.begin(), m2.end());
}

std::ostream& operator<<(std::ostream& os, monomial_t const& m) {
 monomial_t::const_iterator it = std::begin(m), end_it = std::end(m);
 auto next_it = it; ++next_it;
 int power = 1;
 for(;it != end_it; ++it, ++next_it) {
  if(next_it == end_it || *next_it != *it) {
   os << (power > 1 ? "[" : "") << *it << (power > 1 ? "]^" + std::to_string(power) : "");
   power = 1;
  } else
   ++power;
 }
 return os;
}

}
}
