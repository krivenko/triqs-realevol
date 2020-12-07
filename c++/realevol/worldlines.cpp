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

//
// Explicit instantiations
//

template struct worldline_desc_t<1>;
template struct worldline_desc_t<2>;
template struct worldline_desc_t<3>;

template class worldlines_maker<time_expr_operator_t>;
template class worldlines_maker<time_interp_operator_t>;

} // namespace realevol
