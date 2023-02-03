/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013-2023, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include <string>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <variant>

#include <triqs/utility/dressed_iterator.hpp>
#include <triqs/utility/exceptions.hpp>

#include <h5/h5.hpp>

namespace realevol::hilbert_space {

 /// Sequence of indices (`std::vector` of int/string variant objects)
 using indices_t = std::vector<std::variant<int, std::string>>;

 /// Structure of the Green's function
 using gf_struct_t = std::vector<std::pair<std::string, std::vector<std::variant<int, std::string>>>>;

 /// The statistics: Boson or Fermion
 enum statistic_enum {Boson, Fermion};

/// This class represents an ordered set of **indices** of the canonical operators (see [[many_body_operator]]) used to build the Fock states.
/**
 Every element of the set is an arbitrarily long sequence of integers/strings (types can be mixed within one sequence).
 The elements are ordered according to the result of `std::lexicographical_compare`.
 @include triqs/hilbert_space/fundamental_operator_set.hpp
 */
class fundamental_operator_set {
 public:

 /// The set represented as a pair of plain `std::vector` (for fermions and bosons)
 using reduction_t = std::pair<std::vector<indices_t>,std::vector<indices_t>>;

 private:
 using map_t = std::map<indices_t, int>; // the table index <-> n
 map_t fermion_map_index_n;
 map_t boson_map_index_n;

 // internal only
 fundamental_operator_set(std::vector<std::vector<std::string>> const&,
                          std::vector<std::vector<std::string>> const& = {});

 public:

 /// Construct an empty set
 fundamental_operator_set() {}

 /// Construct a set with each stored index being a pair of integers `(i,j)`
 /**
   @param fermion_v Indices of fermions; `i` runs from 0 to `fermion_v.size()-1`; `j` runs from 0 to `fermion_v[i].size()-1` for each `i`
   @param boson_v Indices of bosons; `i` runs from 0 to `boson_v.size()-1`; `j` runs from 0 to `boson_v[i].size()-1` for each `i`
  */
 fundamental_operator_set(std::vector<int> const& fermion_v, std::vector<int> const& boson_v = {}) {
  for (int i = 0; i < fermion_v.size(); ++i)
   for (int j = 0; j < fermion_v[i]; ++j) insert_fermion(i, j);
  for (int i = 0; i < boson_v.size(); ++i)
   for (int j = 0; j < boson_v[i]; ++j) insert_boson(i, j);
 }

 /// Construct from a set of generic index sequences
 /**
   @param s Set of indices for fermions
  */
 template <typename IndexType> fundamental_operator_set(std::set<IndexType> const& fermion_s,
                                                        std::set<IndexType> const& boson_s = {}) {
  for (auto const& i : fermion_s) insert_fermion(i);
  for (auto const& i : boson_s) insert_boson(i);
 }

 /// Construct from two vectors of index sequences
 /**
   @param v Pair of vectors of indices
  */
 explicit fundamental_operator_set(reduction_t const& v) {
  for (auto const& i : v.first) insert_from_indices_t(i, Fermion);
  for (auto const& i : v.second) insert_from_indices_t(i, Boson);
 }

 /// Construct fundamental_operator_set on a pair of GF structure objects
 /**
  @param gf_struct_fermion GF structure object (Fermions)
  @param gf_struct_boson GF structure object (Bosons)
  */
 fundamental_operator_set (gf_struct_t const& gf_struct_fermion, gf_struct_t const& gf_struct_boson) {
  for(auto const& block : gf_struct_fermion)
   for(auto const& inner : block.second)
    insert_fermion(block.first, inner);
  for(auto const& block : gf_struct_boson)
   for(auto const& inner : block.second)
    insert_boson(block.first, inner);
 }

 /// Reduce to a `std::pair<std::vector<indices_t>,std::vector<indices_t>>`
 explicit operator reduction_t() const { return std::make_pair(reverse_map(Fermion),reverse_map(Boson)); }

 /// Insert a new index sequence given as `indices_t`
 /**
   @param ind `indices_t` object
   @param stat Choose between bosonic/fermionic operators
  */
 void insert_from_indices_t(indices_t const& ind, statistic_enum stat = Fermion) {
  map_t & map_index_n = (stat == Fermion ? fermion_map_index_n : boson_map_index_n);
  map_index_n.insert({ind, size(stat)});
  // reorder the indices, which are always given in the order of the indices tuple
  map_t m;
  int i = 0;
  for (auto const& p : map_index_n) m.insert({p.first, i++});
  std::swap(m, map_index_n);
 }

 /// Insert a new index sequence for fermions given as multiple `int`/`std::string` arguments
 template <typename... IndexType> void insert_fermion(IndexType const&... ind) {
  insert_from_indices_t(indices_t{ind...}, Fermion);
 }

 /// Insert a new index sequence for fermions given as multiple `int`/`std::string` arguments
 template <typename... IndexType> void insert_boson(IndexType const&... ind) {
  insert_from_indices_t(indices_t{ind...}, Boson);
 }

 /// Number of elements in this set
 /**
   @return Size of the set
  */
 int size() const { return fermion_map_index_n.size() + boson_map_index_n.size(); }

 /// Number of elements in the boson/fermion subset
 /**
   @param stat Choose between bosonic/fermionic operators
   @return Size of the set
  */
 int size(statistic_enum stat) const {
  return (stat == Fermion ? fermion_map_index_n : boson_map_index_n).size();
 }

 /// Check if a given index sequence is in this set
 /**
   @param t Index sequence to look up
   @param stat Choose between bosonic/fermionic operators
   @return `true` if `t` is in this set
  */
 bool has_indices(indices_t const& t, statistic_enum stat = Fermion) const {
  return (stat == Fermion ? fermion_map_index_n : boson_map_index_n).count(t) == 1;
 }

 /// Request position of a given index sequence for *fermions*
 /**
   @param t Index sequence to look up
   @return Position of the requested index sequence
  */
 int operator[](indices_t const& t) const { return pos(t, Fermion); }

 /// Request position of a given index sequence
 /**
   @param t Index sequence to look up
   @param stat Choose between bosonic/fermionic operators
   @return Position of the requested index sequence
  */
 int pos(indices_t const& t, statistic_enum stat) const {
  try {
   return (stat == Fermion ? fermion_map_index_n : boson_map_index_n).at(t);
  } catch(std::out_of_range &) {
   TRIQS_RUNTIME_ERROR << (stat == Fermion ? "Fermionic" : "Bosonic")
   << " operator with indices (" << t << ") does not belong to this fundamental set!";
  }
 }


 /// Comparison with another fundamental operator set
 bool operator==(fundamental_operator_set const& fops) const {
  return fermion_map_index_n == fops.fermion_map_index_n &&
         boson_map_index_n == fops.boson_map_index_n;
 }

 /// Build and return the reverse map: `int` -> `indices_t`
 /**
   @param stat Choose between bosonic/fermionic operators
   @return The reverse map
  */
 std::vector<indices_t> reverse_map(statistic_enum stat = Fermion) const {
  std::vector<indices_t> r(size(stat));
  map_t const& map_index_n = (stat == Fermion ? fermion_map_index_n : boson_map_index_n);
  for (auto const& x : map_index_n) r[x.second] = x.first;
  return r;
 }

 // Dereference type for const_iterator
 struct _cdress {
  indices_t const& index;
  int linear_index;
  _cdress(typename map_t::const_iterator _it) : index(_it->first), linear_index(_it->second) {}
 };

 /// Constant bidirectional iterator over all stored index sequences. For an iterator `it`, `it->index` gives the `indices_t` object pointed by this iterator, and `it->linear_index` is its position in the set.
 using const_iterator = triqs::utility::dressed_iterator<typename map_t::const_iterator, _cdress>;

 /// Return `const_iterator` to the first element of this set
 /**
   @param stat Choose between bosonic/fermionic operators
   @return Iterator to the first index sequence
  */
 const_iterator begin(statistic_enum stat) const noexcept {
  return (stat == Fermion ? fermion_map_index_n : boson_map_index_n).begin();
 }
 /// Return `const_iterator` to the past-the-end element of this set
 /**
   @param stat Choose between bosonic/fermionic operators
   @return Iterator to the past-the-end element
  */
 const_iterator end(statistic_enum stat) const noexcept {
  return (stat == Fermion ? fermion_map_index_n : boson_map_index_n).end();
 }

 /// Compute a union of two sets
 /**
   @param fops1 First fundamental operator set
   @param fops2 Second fundamental operator set
   @return Union of the sets
  */
 friend fundamental_operator_set merge(fundamental_operator_set const& fops1,
                                       fundamental_operator_set const& fops2) {
  fundamental_operator_set fops(fops1);
  for(auto it = fops2.begin(Fermion); it != fops2.end(Fermion); ++it) fops.insert_fermion(it->index);
  for(auto it = fops2.begin(Boson); it != fops2.end(Boson); ++it) fops.insert_boson(it->index);
  return fops;
 }

 /// Write this set as an HDF5 attribute
 /**
   @param obj The HDF5 object the attribute is attached to
   @param name Name of the attribute
   @param f Fundamental set to write
  */

 friend void h5_write_attribute(h5::object obj, std::string const& name, fundamental_operator_set const& f);
 /// Read a set from an HDF5 attribute
 /**
   @param obj The HDF5 object the attribute is attached to
   @param name Name of the attribute
   @param f Reference to a fundamental set to be read
  */
 friend void h5_read_attribute(h5::object obj, std::string const& name, fundamental_operator_set& f);

};
}
