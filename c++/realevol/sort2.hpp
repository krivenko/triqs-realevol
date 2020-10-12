/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 T. Shields
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

// Sort two vectors in the same way, with criteria that uses only one of the vectors.
// Implementation based on
// http://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of

#include <vector>
#include <algorithm>
#include <numeric>

namespace realevol {

// Find a sort permutation
template<typename V, typename Compare = std::less<typename V::value_type>>
std::vector<std::size_t> sort_permutation(const V& vec, Compare comp = {}) {
 std::vector<std::size_t> p(vec.size());
 std::iota(p.begin(), p.end(), 0);
 std::sort(p.begin(), p.end(),
  [&](std::size_t i, std::size_t j){ return comp(vec[i], vec[j]); }
 );
 return p;
}

// Apply a sort permutation
template<typename V>
V apply_permutation(V const& vec, std::vector<std::size_t> const& p) {
 V sorted(vec.size());
 std::transform(p.begin(), p.end(), sorted.begin(),
                [&](std::size_t i){ return vec[i]; }
 );
 return sorted;
}

template<typename V>
void apply_permutation_in_place(V & vec, std::vector<std::size_t> const& p) {
 std::vector<bool> done(vec.size());
 for (std::size_t i = 0; i < vec.size(); ++i) {
  if (done[i]) continue;
  done[i] = true;
  std::size_t prev_j = i;
  std::size_t j = p[i];
  while (i != j) {
   using std::swap;
   swap(vec[prev_j], vec[j]);
   done[j] = true;
   prev_j = j;
   j = p[j];
  }
 }
}

}
