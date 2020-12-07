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
#include <triqs/utility/first_include.hpp>

#include <string>

#include "init_state.hpp"

namespace realevol {

std::ostream & operator<<(std::ostream & os, init_state const& st) {
 os << "Fundamental operator set:" << std::endl;
 os << " Fermions: ";
 for(auto it = st.fops.begin(Fermion); it != st.fops.end(Fermion); ++it)
  os << "(" << it->index << ") ";
 os << std::endl << " Bosons: ";
 for(auto it = st.fops.begin(Boson); it != st.fops.end(Boson); ++it)
  os << "(" << it->index << ") ";

 os << std::endl << "Dimension of fermionic space: " << st.full_hs.size(Fermion) << std::endl;
 os << "Dimension of bosonic space: " << st.full_hs.size(Boson) << std::endl;

 os << st.sub_hilbert_spaces.size() << " relevant invariant subspaces:" << std::endl;
 for(auto const& sp : st.sub_hilbert_spaces) {
  os << " Subspace " << sp.get_index() << " [" << sp.size() << "]: ";
  for(auto f : sp.get_all_fock_states()) os << f << " ";
  os << std::endl;
 }
 os << std::endl;

 for(int i = 0; i < st.weighted_states.size(); ++i) {
  auto const& s = st.weighted_states[i];
  os << "State " << i << " within subspace " << s.state.get_hilbert().get_index()
     << ", weight = " << s.weight << ":" << s.state << std::endl;
 }

 return os;
}

// -----------------------------------------------------------------

void h5_write(h5::group fg, std::string const &name, init_state const& st) {
 auto gr = fg.create_group(name);
 write_hdf5_format(gr, st);

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

  double weight = {};
  h5_read(wst_gr, "weight", weight);
  int sp_index = {};
  h5_read(wst_gr, "sp_index", sp_index);
  st.weighted_states.emplace_back(state_on_subspace_t(st.sub_hilbert_spaces[sp_index]), weight);
  h5_read(wst_gr, "amplitudes", st.weighted_states.back().state.amplitudes());
 }
}

}
