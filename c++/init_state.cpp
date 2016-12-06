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
#include <triqs/utility/first_include.hpp>

#include <cmath>
#include <triqs/arrays/vector.hpp>
#include <arpack/arpack_worker.hpp>

#include "init_state.hpp"

namespace realevol {

std::ostream & operator<<(std::ostream & os, init_state const& st) {
 os << "Fundamental operator set" << std::endl;
 os << "------------------------" << std::endl;
 os << "Fermions: ";
 for(auto it = st.fops.begin(Fermion); it != st.fops.end(Fermion); ++it)
  os << "(" << it->index << ") ";
 os << std::endl << "Bosons: ";
 for(auto it = st.fops.begin(Boson); it != st.fops.end(Boson); ++it)
  os << "(" << it->index << ") ";

 os << std::endl << "Dimension of fermionic space: " << st.full_hs.size(Fermion) << std::endl;
 os << "Dimension of bosonic space: " << st.full_hs.size(Boson) << std::endl;

 os << st.sub_hilbert_spaces.size() << " relevant invariant subspaces with dimensions ";
 for(auto const& sp : st.sub_hilbert_spaces) os << sp.size() << " ";
 os << std::endl;

 for(int i = 0; i < st.weighted_states.size(); ++i) {
  auto const& s = st.weighted_states[i];
  os << "State " << i << " within subspace " << s.state.get_hilbert().get_index()
     << ", weight = " << s.weight << ":" << std::endl << s.state;
 }

 return os;
}

// -----------------------------------------------------------------

void h5_write(h5::group fg, std::string const &name, init_state const& st) {
 auto gr = fg.create_group(name);
 gr.write_triqs_hdf5_data_scheme(st);

 h5_write_attribute(gr, "fops", st.fops);
 h5_write(gr, "full_hs", st.full_hs);
 h5_write(gr, "sub_hilbert_spaces", st.sub_hilbert_spaces);

 auto states_gr = gr.create_group("weighted_states");
 for(int i = 0; i < st.weighted_states.size(); ++i) {
  auto wst_gr = states_gr.create_group(std::to_string(i));
  auto const& wst = st.weighted_states[i];

  h5_write(wst_gr, "weight", wst.weight);
  h5_write(wst_gr, "sp_index", wst.state.get_hilbert().get_index());
  h5_write(wst_gr, "amplitudes", wst.state.amplitudes());
 }
}

// -----------------------------------------------------------------

void h5_read(h5::group fg, std::string const &name, init_state & st) {
 auto gr = fg.open_group(name);

 h5_read_attribute(gr, "fops", st.fops);
 h5_read(gr, "full_hs", st.full_hs);
 h5_read(gr, "sub_hilbert_spaces", st.sub_hilbert_spaces);

 auto states_gr = gr.open_group("weighted_states");
 int n_wst = states_gr.get_all_subgroup_names().size();
 for(int i = 0; i < n_wst; ++i) {
  auto wst_gr = states_gr.open_group(std::to_string(i));

  double weight;
  h5_read(wst_gr, "weight", weight);
  unsigned long sp_index;
  h5_read(wst_gr, "sp_index", sp_index);
  st.weighted_states.emplace_back(state_on_subspace_t(st.sub_hilbert_spaces[sp_index]), weight);
  h5_read(wst_gr, "amplitudes", st.weighted_states.back().state.amplitudes());
 }
}

// -----------------------------------------------------------------

static_operator_t make_static_op(operator_t const& op, std::string const& error_message) {
 static_operator_t static_op;
 for(auto const& m : op) {
  if(!is_constant(m.coef)) TRIQS_RUNTIME_ERROR << error_message;
  auto new_coef = m.coef(0);
  static_operator_t new_monomial(new_coef);
  for(auto const& c : m.monomial)
   new_monomial *= static_operator_t::make_canonical(c.stat, c.dagger, c.indices);
  static_op += new_monomial;
 }
 return static_op;
}

init_state make_pure_init_state(operator_t const& generator,
                                fundamental_operator_set const& fops,
                                std::map<operators::indices_t, int> const& bits_per_boson) {

 // Static version of generator
 auto gen = make_static_op(generator, "Generating operator must be time-independent!");

 init_state st(fops, bits_per_boson);

 // Add one subspace with index 0
 st.sub_hilbert_spaces.emplace_back(0);
 auto & sp = st.sub_hilbert_spaces[0];

 // Vaccum state in the full Hilbert space
 state_on_space_t vacuum(st.full_hs);
 vacuum(fock_state_t(0)) = 1;

 // Imperative version of generator
 static_op_on_space_t imp_gen(gen, st.fops, st.full_hs);

 auto psi = imp_gen(vacuum);

 // Fill st.sub_hilbert_spaces[0]
 foreach(psi, [&sp](fock_state_t f, dcomplex){ sp.add_fock_state(f); });

 // Add a weighted state
 state_on_subspace_t st_on_sp(sp);
 foreach(psi, [&st_on_sp,&sp](fock_state_t f, dcomplex a){ st_on_sp(sp.get_state_index(f)) = a; });
 st.weighted_states.emplace_back(std::move(st_on_sp), 1.0);

 return st;
}

/*
std::vector<init_state> init_statehermal(operator_t const& h0, hilbert_space_structure & hss,
                                             double beta, double temperature_cutoff) {

 // Static version the equilibrium Hamiltonian is static
 auto h0_ = make_static_op(h0, "Equilibrium Hamiltonian must be time-independent!");

 // Hilbert space structure generated by the equilibrium Hamiltonian
 hilbert_space_structure hss0(h0, hss.fops, hss.full_hs, false);

 // Time complexity on Lanczos iterations: O(k(nnz + n))
 // n - matrix size, nnz - number of non-zero elements, k - number of Krylov basis vectors
 // J. Chen and Y. Saad, "Lanczos Vectors versus Singular Vectors for Effective Dimension Reduction,"
 // IEEE Transactions on Knowledge and Data Engineering, vol. 21, no. 8, pp. 1091-1103, Aug. 2009.

 // Iterate over all sectors of h0
//  for(auto const& sp : hss0.sub_hilbert_spaces) {
//   static_op_on_subspace_t imp_h0(h0_, hss0.fops, hss.full_hs);
//
//   auto lw = lanczos_worker<static_op_on_subspace_t, state_on_subspace_t>(imp_h0,1e-10);
//
//   state_on_subspace_t psi0(sp); psi0(0) = 1.0;
//   lw(psi0);
//
//   std::cout << lw.values() << std::endl;
//  }

 std::vector<init_state> states;
 // TODO

 return states;
}
*/
}
