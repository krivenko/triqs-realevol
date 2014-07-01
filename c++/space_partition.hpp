#pragma once

#include <set>
#include <map>
#include <utility>
#include <limits>
#include <boost/pending/disjoint_sets.hpp>

namespace realevol {

// Requirements
// ------------
// * State must model StateVector
// * StateType OperatorType::operator()(StateType const &st) must apply the operator to 'st' and return the result

template <typename StateType, typename OperatorType> class space_partition {

 public:
 using index_t = uint32_t;
 using state_t = StateType;
 using operator_t = OperatorType;
 using amplitude_t = typename state_t::value_type;

 using block_mapping_t = std::set<std::pair<index_t,index_t>>;

 static constexpr amplitude_t tolerance = std::numeric_limits<amplitude_t>::epsilon();

 space_partition(state_t const& st, operator_t const& H)
    : subspaces(st.size()), tmp_state(make_zero_state(st)) {

  // Iteration over all initial basis states
  for (index_t i = 0; i < tmp_state.size(); ++i) {
   state_t initial_state = tmp_state;
   initial_state(i) = amplitude_t(1);

   state_t final_state = H(initial_state);

   // Iterate over non-zero final amplitudes
   foreach(final_state, [&](index_t f, amplitude_t amplitude) {
    if (std::fabs(amplitude) < tolerance) return;
    auto i_subspace = subspaces.find_set(i);
    auto f_subspace = subspaces.find_set(f);
    if (i_subspace != f_subspace) subspaces.link(i_subspace, f_subspace);
   });
  }

  _update_index();
 }

 space_partition(space_partition const&) = default;

 // Access information about subspaces
 index_t n_subspaces() const { return representative_to_index.size(); }

 template <typename Lambda> friend void foreach(space_partition& SP, Lambda L) {
  for (index_t n = 0; n < SP.tmp_state.size(); ++n) L(n, SP.lookup_basis_state(n));
 };

 index_t lookup_basis_state(index_t basis_state) { return representative_to_index[subspaces.find_set(basis_state)]; }

 block_mapping_t find_mappings(operator_t const& op, bool diagonal_only = false) {

  block_mapping_t mapping;

  // Iteration over all initial basis states
  for (index_t i = 0; i < tmp_state.size(); ++i) {
   state_t initial_state = tmp_state;
   initial_state(i) = amplitude_t(1);
   auto i_subspace = subspaces.find_set(i);

   state_t final_state = op(initial_state);

   // Iterate over non-zero final amplitudes
   foreach(final_state, [&](index_t f, amplitude_t amplitude) {
    if (std::fabs(amplitude) < tolerance) return;
    auto f_subspace = subspaces.find_set(f);
    if((!diagonal_only) || i_subspace==f_subspace)
      mapping.insert(std::make_pair(i_subspace,f_subspace));
   });
  }

  return mapping;
 }

 private:
 void _update_index() {
  auto p = subspaces.parents();
  subspaces.compress_sets(p.begin(), p.end());  // parents are representatives
  subspaces.normalize_sets(p.begin(), p.end()); // the representative has the smallest index in the set

  // Update representative_to_index
  representative_to_index.clear();
  for (index_t n = 0; n < tmp_state.size(); ++n) {
   representative_to_index.insert(std::make_pair(subspaces.find_set(n), representative_to_index.size()));
  }
 }

 // Temporary zero state
 mutable state_t tmp_state;
 // Subspaces
 boost::disjoint_sets_with_storage<> subspaces;
 // Map representative basis state to subspace index
 std::map<index_t, index_t> representative_to_index;
};
}
