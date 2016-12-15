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

#include "common.hpp"

#include <vector>
#include <set>
#include <triqs/hilbert_space/space_partition.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/gfs.hpp>

namespace realevol {

using triqs::arrays::matrix;

struct hilbert_space_structure {

 // Fundamental operator set
 fundamental_operator_set const& fops;
 // Full hilbert space
 class hilbert_space const& full_hs;
 // Invariant subspaces
 std::vector<sub_hilbert_space> sub_hilbert_spaces;
 // Connections between subspaces: *_connection[stat][operator_linear_index][B] -> B'
 matrix<long> creation_connection[2], annihilation_connection[2];
 // Hilbert space partition
 space_partition<dyn_state_on_space_t, op_on_space_t> partition;

 /// Partition a space
 template<typename IsZeroPredicate>
 hilbert_space_structure(operator_t const& h,
                         fundamental_operator_set const& fops,
                         class hilbert_space const& full_hs,
                         fundamental_operator_set const& merger_fops,
                         IsZeroPredicate is_zero_pred) :
  fops(fops), full_hs(full_hs),
  partition(dyn_state_on_space_t(full_hs), op_on_space_t(h, fops, full_hs), false, is_zero_pred) {
  merge(merger_fops, is_zero_pred);
  fill_subspaces();
 }

 /// For a given subspace sp, returns a set of indices of sub_hilbert_spaces, sp has basis states in.
 std::set<long> subspace_branching(sub_hilbert_space const& sp) {
  std::set<long> branching;
  for(fock_state_t f : sp.get_all_fock_states()) {
   branching.insert(partition.lookup_basis_state(f));
  }
  return branching;
 }

private:

 template<typename IsZeroPredicate>
 void merge(fundamental_operator_set const& merger_fops, IsZeroPredicate is_zero_pred) {

  using realevol::hilbert_space::Fermion;
  using realevol::hilbert_space::Boson;

  // Merge subspaces
  std::vector<op_on_space_t> create_ops[2], destroy_ops[2];
  for(auto stat : {Fermion, Boson}) {
   create_ops[stat].reserve(merger_fops.size(stat));
   destroy_ops[stat].reserve(merger_fops.size(stat));

   for (auto it = merger_fops.begin(stat); it != merger_fops.end(stat); ++it) {
    auto create = operator_t::make_canonical(stat, true, it->index);
    auto destroy = operator_t::make_canonical(stat, false, it->index);

    create_ops[stat].emplace_back(create, merger_fops, full_hs);
    destroy_ops[stat].emplace_back(destroy, merger_fops, full_hs);

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
    int n = it->linear_index;

    auto create_conns = partition.find_mappings(create_ops[stat][n], false, is_zero_pred);
    for(auto const& conn : create_conns) creation_connection[stat](n, conn.first) = conn.second;

    auto destroy_conns = partition.find_mappings(destroy_ops[stat][n], false, is_zero_pred);
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
