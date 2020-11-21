/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <set>
#include <vector>

#include <triqs/arrays/matrix.hpp>
#include <triqs/gfs.hpp>

#include "hilbert_space/space_partition.hpp"

#include "common.hpp"

namespace realevol {

using triqs::arrays::matrix;

/// List of subspace branchings
///      /-> sp'_1              /-> sp'_1
/// sp_1 --> sp'_2 , ... , sp_N --> sp'_2
///      \-> sp'_3              \-> sp'_3
using branchings_t = std::vector<std::set<long>>;

template<typename HamiltonianType>
struct hilbert_space_structure {

 using HScalarType = typename HamiltonianType::scalar_t;

 // Fundamental operator set
 fundamental_operator_set const& fops;
 // Full hilbert space
 class hilbert_space const& full_hs;
 // Invariant subspaces
 std::vector<sub_hilbert_space> sub_hilbert_spaces;
 // Connections between subspaces: *_connection[stat][operator_linear_index][B] -> B'
 // WARNING: For operators belonging to fops but missing from merger_fops
 // (not 1-to-1 mappings) no connections are stored
 matrix<long> creation_connection[2], annihilation_connection[2];
 // Hilbert space partition
 mutable space_partition<dyn_state_on_space_t<HScalarType>, op_on_space_t<HScalarType>> partition;

 /// Partition a space
 template<typename IsZeroPredicate>
 hilbert_space_structure(HamiltonianType const& h,
                         fundamental_operator_set const& fops,
                         class hilbert_space const& full_hs,
                         fundamental_operator_set const& merger_fops,
                         IsZeroPredicate is_zero_pred) :
  fops(fops), full_hs(full_hs),
  partition(dyn_state_on_space_t<HScalarType>(full_hs),
            op_on_space_t<HScalarType>(h, fops, full_hs),
            false,
            is_zero_pred) {
  merge(merger_fops, is_zero_pred);
  fill_subspaces();
 }

 /// For a given subspace sp, returns a set of indices of sub_hilbert_spaces, sp has basis states in.
 std::set<long> compute_branching(sub_hilbert_space const& sp) const {
  std::set<long> branching;
  for(fock_state_t f : sp.get_all_fock_states()) {
   branching.insert(partition.lookup_basis_state(f));
  }
  return branching;
 }

 /// For a given subspace sp, returns a set of indices of sub_hilbert_spaces, sp has basis states in.
 branchings_t compute_branchings(std::vector<sub_hilbert_space> const& subspaces) const {
  branchings_t branchings;
  branchings.reserve(subspaces.size());
  for(auto const& sp : subspaces) branchings.emplace_back(compute_branching(sp));
  return branchings;
 }

 template<typename ClassifyPredicate>
 std::vector<bool> classify_subspaces(HamiltonianType const& op, ClassifyPredicate classify_pred) const {
  std::vector<bool> res(sub_hilbert_spaces.size());

  dyn_state_on_space_t<HScalarType> from_st(full_hs), to_st(full_hs);
  op_on_space_t<HScalarType> imp_op(op, fops, full_hs);

  for(long spn = 0; spn < sub_hilbert_spaces.size(); ++spn) {
   auto const& hs = sub_hilbert_spaces[spn];
   res[spn] = true;
   for(auto f : hs.get_all_fock_states()) {
    from_st(f) = 1;
    imp_op.apply(from_st, to_st);
    from_st.zero();
    if(!classify_pred(to_st)) {
     res[spn] = false;
     break;
    }
   }
  }
  return res;
 }

private:

 template<typename IsZeroPredicate>
 void merge(fundamental_operator_set const& merger_fops, IsZeroPredicate is_zero_pred) {

  using realevol::hilbert_space::Fermion;
  using realevol::hilbert_space::Boson;

  // Merge subspaces
  std::vector<op_on_space_t<HScalarType>> create_ops[2], destroy_ops[2];
  for(auto stat : {Fermion, Boson}) {
   create_ops[stat].reserve(merger_fops.size(stat));
   destroy_ops[stat].reserve(merger_fops.size(stat));

   for (auto it = merger_fops.begin(stat); it != merger_fops.end(stat); ++it) {
    auto create = HamiltonianType::make_canonical(stat, true, it->index);
    auto destroy = HamiltonianType::make_canonical(stat, false, it->index);

    create_ops[stat].emplace_back(create, fops, full_hs);
    destroy_ops[stat].emplace_back(destroy, fops, full_hs);

    partition.merge_subspaces(create_ops[stat].back(), destroy_ops[stat].back(), true, is_zero_pred);
   }
  }

  // Fill connections, only for operators in merger_fops
  // (other operators are not necessarily 1-to-1 mappings)
  for(auto stat : {Fermion, Boson}) {
   creation_connection[stat].resize(fops.size(stat), partition.n_subspaces());
   annihilation_connection[stat].resize(fops.size(stat), partition.n_subspaces());

   creation_connection[stat].as_array_view() = -1;
   annihilation_connection[stat].as_array_view() = -1;

   for (auto it = merger_fops.begin(stat); it != merger_fops.end(stat); ++it) {
    int n = fops.pos(it->index, stat);

    auto create_conns = partition.find_mappings(create_ops[stat][it->linear_index], false, is_zero_pred);
    for(auto const& conn : create_conns) creation_connection[stat](n, conn.first) = conn.second;

    auto destroy_conns = partition.find_mappings(destroy_ops[stat][it->linear_index], false, is_zero_pred);
    for(auto const& conn : destroy_conns) annihilation_connection[stat](n, conn.first) = conn.second;
   }
  }
 }

 // Fill subspaces
 void fill_subspaces() {
  sub_hilbert_spaces.reserve(partition.n_subspaces());
  for (int n = 0; n < partition.n_subspaces(); ++n) sub_hilbert_spaces.emplace_back(n);
  foreach(partition, [&](fock_state_t s, int spn) { sub_hilbert_spaces[spn].add_fock_state(s); });
 }

};

}
