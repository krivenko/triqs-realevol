/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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
#include "worldlines.hpp"

namespace realevol {

//
// class worldlines_maker
//

template<typename HamiltonianType>
std::vector<worldline_desc_t<2>>
worldlines_maker<HamiltonianType>::make_gf_worldlines(
  gf_struct_t const& gf_struct,
  bool is_greater) {
  auto const& c_conn = hss.annihilation_connection[Fermion];
  auto const& cdag_conn = hss.creation_connection[Fermion];

  std::vector<worldline_desc_t<2>> res;
  auto const& fops = initial_state.get_fops();
  auto const& wst = initial_state.get_weighted_states();

  int block_index = 0;
  for(auto const& block : gf_struct) {
    int bl_size = block.second.size();
    for (int inner_index1 = 0; inner_index1 < bl_size; ++inner_index1) {
      indices_t index1 = {block.first, block.second[inner_index1]};
      auto M_c = monomial_t{{Fermion, false, index1}};
      int n1 = fops[index1]; // linear_index of c
      for (int inner_index2 = 0; inner_index2 < bl_size; ++inner_index2) {
        indices_t index2 = {block.first, block.second[inner_index2]};
        auto M_c_dag = monomial_t{{Fermion, true, index2}};
        int n2 = fops[index2]; // linear_index of c^+

        for(int wst_i = 0; wst_i < wst.size(); ++wst_i) {
          auto const& st = wst[wst_i];
          auto const& branching = branchings[st.state.get_hilbert().get_index()];
          for(long right_sp_index : branching) {
            // Checking selection rule for <middle_sp_index| |right_sp_index>
            long middle_sp_index = is_greater ? cdag_conn(n2,right_sp_index) :
                                                c_conn(n1,right_sp_index);
            if(middle_sp_index == -1) continue;
            // Checking selection rule for <left_sp_index| |middle_sp_index>
            long left_sp_index = is_greater ? c_conn(n1,middle_sp_index) :
                                              cdag_conn(n2,middle_sp_index);
            if(left_sp_index == -1 || !branching.count(left_sp_index)) continue;

            if(is_greater) {
              //  Greater component G^>_{index1,index2}(t,t'):
              //   (-i/hbar) * <l| c_index1(t) |m><m| c^+_index2(t') |r> * weight
              res.push_back(
                worldline_desc_t<2>{
                  worldline_desc_t<2>::GreaterGf, block_index, inner_index1, inner_index2,
                  -1i / hbar,
                  std::array<monomial_t, 2>{M_c_dag, M_c},
                  std::array<long, 3>{right_sp_index, middle_sp_index, left_sp_index},
                  wst_i
                }
              );
            } else {
              //  Lesser component G^<_{index1,index2}(t,t'):
              //   (i/hbar) * <l| c^+_index2(t') |m><m| c_index1(t) |r> * weight
              res.push_back(
                worldline_desc_t<2>{
                  worldline_desc_t<2>::LesserGf, block_index, inner_index1, inner_index2,
                  1i / hbar,
                  std::array<monomial_t, 2>{M_c, M_c_dag},
                  std::array<long, 3>{right_sp_index, middle_sp_index, left_sp_index},
                  wst_i
                }
              );
            }
          }
        }
      }
    }
    ++block_index;
  }

  return res;
}

template<typename HamiltonianType>
std::vector<worldline_desc_t<2>>
worldlines_maker<HamiltonianType>::make_chi_worldlines(
  chi_indices_t const& chi_indices
) {
  auto const& c_conn = hss.annihilation_connection[Fermion];

  std::vector<worldline_desc_t<2>> res;
  auto const& fops = initial_state.get_fops();
  auto const& wst = initial_state.get_weighted_states();

  int chi_size = chi_indices.size();
  for (int inner_index1 = 0; inner_index1 < chi_size; ++inner_index1) {
    indices_t index1 = {chi_indices[inner_index1].first,
                        chi_indices[inner_index1].second};
    auto M_n_1 = monomial_t{{Fermion, true, index1}, {Fermion, false, index1}};
    int n1 = fops[index1]; // linear_index of n1
    for (int inner_index2 = 0; inner_index2 < chi_size; ++inner_index2) {
      indices_t index2 = {chi_indices[inner_index2].first,
                          chi_indices[inner_index2].second};
      auto M_n_2 = monomial_t{{Fermion, true, index2}, {Fermion, false, index2}};
      int n2 = fops[index2]; // linear_index of n2

      for(int wst_i = 0; wst_i < wst.size(); ++wst_i) {
        auto const& st = wst[wst_i];
        auto const& branching = branchings[st.state.get_hilbert().get_index()];
        for(long right_sp_index : branching) {
          // Checking selection rule for <middle_sp_index| |right_sp_index>
          long middle_sp_index = (c_conn(n2,right_sp_index) != -1) ? right_sp_index : -1;
          if(middle_sp_index == -1) continue;
          // Checking selection rule for <left_sp_index| |middle_sp_index>
          long left_sp_index = (c_conn(n1,middle_sp_index) != -1) ? middle_sp_index : -1;
          if(left_sp_index == -1 || !branching.count(left_sp_index)) continue;

          //  Susceptibility \chi_{index1,index2}(t,t'):
          //   (-i/hbar) * <l| n_index1(t) |m><m| n_index2(t') |r> * weight
          res.push_back(
            worldline_desc_t<2>{
              worldline_desc_t<2>::Susceptibility, 0, inner_index1, inner_index2,
              -1i / hbar,
              std::array<monomial_t, 2>{M_n_2, M_n_1},
              std::array<long, 3>{right_sp_index, middle_sp_index, left_sp_index},
              wst_i
            }
          );
        }
      }
    }
  }

  return res;
}

template<typename HamiltonianType>
template<std::size_t NPoints, std::size_t Point>
void worldlines_maker<HamiltonianType>::make_worldlines_impl(
  std::array<std::pair<op_iter_t, op_iter_t>, NPoints> op_iters,
  std::vector<worldline_desc_t<NPoints>> & worldlines
) const {
  if constexpr(Point + 1 < NPoints) {
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
          worldline_desc_t<NPoints>::GreaterGf, // TODO: Remove
          0,                                    // TODO: Remove
          0,                                    // TODO: Remove
          0,                                    // TODO: Remove
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

template struct worldline_desc_t<1>;
template struct worldline_desc_t<2>;
template struct worldline_desc_t<3>;

template class worldlines_maker<time_expr_operator_t>;    // TODO: Remove
template class worldlines_maker<time_interp_operator_t>;  // TODO: Remove

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
