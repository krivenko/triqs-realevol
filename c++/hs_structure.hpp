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

#include <triqs/arrays/matrix.hpp>
#include "triqs/operators/many_body_operator.hpp"
#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include <triqs/hilbert_space/space_partition.hpp>
#include <triqs/gfs.hpp>

#include "time_expr.hpp"

namespace realevol {

using triqs::arrays::matrix;
using triqs::utility::is_zero;
namespace hsns = realevol::hilbert_space; // FIXME

struct is_zero_on_mesh {
 triqs::gfs::segment_mesh const& mesh;
 is_zero_on_mesh(triqs::gfs::segment_mesh const& mesh) : mesh(mesh) {}
 bool operator()(time_expr const& te) {
  if(te.is_zero()) return true;
  for(auto t : mesh)
   if(!triqs::utility::is_zero(te(t))) return false;
  return true;
 }
};

struct hilbert_space_structure {

 using operator_t = realevol::operators::many_body_operator_generic<time_expr>;

 // Full hilbert space
 hsns::hilbert_space full_hs;
 // Invariant subspaces
 std::vector<hsns::sub_hilbert_space> sub_hilbert_spaces;
 // Connections between subspaces: *_connection[stat][operator_linear_index][B] -> B'
 matrix<long> creation_connection[2], annihilation_connection[2];

 template<typename IsZeroPredicate>
 hilbert_space_structure(hsns::fundamental_operator_set const& fops,
                         operator_t const& h, std::vector<int> bits_per_boson,
                         bool merge_subspaces,
                         IsZeroPredicate is_zero_pred) :
  full_hs(fops, bits_per_boson) {

  using hsns::hilbert_space;
  using hsns::fock_state_t;
  using hsns::state;
  using hsns::imperative_operator;
  using hsns::space_partition;

  using imp_operator_t = imperative_operator<hilbert_space, time_expr, false>;

  imp_operator_t hamiltonian(h, fops, full_hs);
  state<hilbert_space, time_expr, true> st(full_hs);

  using space_partition_t = space_partition<state<hilbert_space, time_expr, true>,
                                            imp_operator_t, IsZeroPredicate>;
  using realevol::hilbert_space::Fermion;
  using realevol::hilbert_space::Boson;

  // Split the Hilbert space
  space_partition_t SP(st, hamiltonian, false, is_zero_pred);

  // Merge subspaces
  std::vector<imp_operator_t> create_ops[2], destroy_ops[2];
  if(merge_subspaces) {
   for(auto stat : {Fermion, Boson}) {
    create_ops[stat].reserve(fops.size(stat));
    destroy_ops[stat].reserve(fops.size(stat));

    for (auto it = fops.begin(stat); it != fops.end(stat); ++it) {
     auto create = operator_t::make_canonical(stat, true, it->index);
     auto destroy = operator_t::make_canonical(stat, false, it->index);

     create_ops[stat].emplace_back(create, fops, full_hs);
     destroy_ops[stat].emplace_back(destroy, fops,  full_hs);

     SP.merge_subspaces(create_ops[stat].back(), destroy_ops[stat].back(), true);
    }
   }
  }

  // Fill subspaces
  sub_hilbert_spaces.reserve(SP.n_subspaces());
  for (int n = 0; n < SP.n_subspaces(); ++n) sub_hilbert_spaces.emplace_back(n);

  foreach(SP, [&](fock_state_t s, int spn) { sub_hilbert_spaces[spn].add_fock_state(s); });

  if(!merge_subspaces) return;

  // Fill connections
  for(auto stat : {Fermion, Boson}) {
   creation_connection[stat].resize(fops.size(stat), SP.n_subspaces());
   annihilation_connection[stat].resize(fops.size(stat), SP.n_subspaces());

   creation_connection[stat].as_array_view() = -1;
   annihilation_connection[stat].as_array_view() = -1;

   for (auto it = fops.begin(stat); it != fops.end(stat); ++it) {
    int n = it->linear_index;

    auto create_conns = SP.find_mappings(create_ops[stat][n], false);
    for(auto const& conn : create_conns) creation_connection[stat](n, conn.first) = conn.second;

    auto destroy_conns = SP.find_mappings(destroy_ops[stat][n], false);
    for(auto const& conn : destroy_conns) annihilation_connection[stat](n, conn.first) = conn.second;
   }
  }
 }

};

}
