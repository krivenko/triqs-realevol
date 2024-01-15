/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2024, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <numeric>
#include <utility>

#include "array_utility.hpp"
#include "time_expr.hpp"
#include "time_interp.hpp"
#include "worldlines.hpp"

namespace realevol {

template<typename HamiltonianType>
template<std::size_t NPoints, std::size_t Point>
void worldlines_maker<HamiltonianType>::make_worldlines_impl(
  std::array<std::pair<op_iter_t, op_iter_t>, NPoints> op_iters,
  std::vector<worldline_desc_t<NPoints>> & worldlines
) const {
  if constexpr(Point < NPoints) {
    auto & [it, end_it] = op_iters[Point];
    for(; it != end_it; ++it) { // Iterate over monomials in ops[Point]
      make_worldlines_impl<NPoints, Point + 1>(op_iters, worldlines);
    }
  } else {
    dcomplex coef = std::accumulate(
      op_iters.begin(),
      op_iters.end(),
      dcomplex(1.0),
      [](dcomplex prod, std::pair<op_iter_t, op_iter_t> const& it) {
        return prod * dcomplex(it.first->coef);
      }
    );

    auto M = map_array<monomial_t>(
      [](std::pair<op_iter_t, op_iter_t> const& it) {
        return it.first->monomial;
      },
      op_iters
    );

    auto const& wst = initial_state.get_weighted_states();

    auto make_sp_indices = [&](long right_sp_index, std::set<long> const& branching) ->
      std::optional<std::array<long, NPoints + 1>> {
      std::array<long, NPoints + 1> sp_indices;
      sp_indices[0] = right_sp_index;
      for(int p = 0; p < NPoints; ++p) {
        sp_indices[p + 1] = hss.make_monomial_connections(M[p])[sp_indices[p]];
        if(sp_indices[p + 1] == -1) return {};
      }
      if(!branching.count(sp_indices[NPoints])) return {};
      return sp_indices;
    };

    for(int wst_i = 0; wst_i < wst.size(); ++wst_i) {
      auto const& st = wst[wst_i];
      auto const& branching = branchings[st.state.get_hilbert().get_index()];

      for(long right_sp_index : branching) {

        auto sp_indices = make_sp_indices(right_sp_index, branching);
        if(!sp_indices) continue;

        worldlines.emplace_back(worldline_desc_t<NPoints>{
          coef,
          M,
          std::move(sp_indices.value()),
          wst_i
        });
      }
    }
  }
}

template<typename HamiltonianType>
template<std::size_t NPoints>
std::vector<worldline_desc_t<NPoints>>
worldlines_maker<HamiltonianType>::make_worldlines(
  std::array<static_operator_t, NPoints> const& ops
) const {
  std::vector<worldline_desc_t<NPoints>> worldlines;
  auto op_iters = map_array<std::pair<op_iter_t, op_iter_t>>(
    [](auto const& op) { return std::make_pair(op.begin(), op.end()); },
  ops);
  make_worldlines_impl<NPoints, 0>(std::move(op_iters), worldlines);
  return worldlines;
}

//
// Explicit instantiations
//

using time_expr_operator_t = realevol::operators::many_body_operator_generic<time_expr>;
using time_interp_operator_t = realevol::operators::many_body_operator_generic<time_interp>;

template struct worldline_desc_t<1>;
template struct worldline_desc_t<2>;
template struct worldline_desc_t<3>;

template std::vector<worldline_desc_t<1>>
worldlines_maker<time_expr_operator_t>::make_worldlines<1>(std::array<static_operator_t, 1> const&) const;
template std::vector<worldline_desc_t<2>>
worldlines_maker<time_expr_operator_t>::make_worldlines<2>(std::array<static_operator_t, 2> const&) const;
template std::vector<worldline_desc_t<3>>
worldlines_maker<time_expr_operator_t>::make_worldlines<3>(std::array<static_operator_t, 3> const&) const;

template std::vector<worldline_desc_t<1>>
worldlines_maker<time_interp_operator_t>::make_worldlines<1>(std::array<static_operator_t, 1> const&) const;
template std::vector<worldline_desc_t<2>>
worldlines_maker<time_interp_operator_t>::make_worldlines<2>(std::array<static_operator_t, 2> const&) const;
template std::vector<worldline_desc_t<3>>
worldlines_maker<time_interp_operator_t>::make_worldlines<3>(std::array<static_operator_t, 3> const&) const;

} // namespace realevol
