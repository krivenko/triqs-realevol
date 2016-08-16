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

#include <algorithm>

#include "solver.hpp"
#include "hs_structure.hpp"
#include "init_state.hpp"

namespace realevol {

struct index_visitor  {
 std::vector<std::string> indices;
 void operator()(int i) { indices.push_back(std::to_string(i)); }
 void operator()(std::string s) { indices.push_back(s); }
};

solver::solver(std::map<std::string,indices_type> const& gf_struct, std::pair<double,double> time_window, int n_t) :
  gf_struct(gf_struct) {

 std::vector<std::string> block_names;
 std::vector<g_2t_t> g_ret_blocks;
 std::vector<g_2t_t> g_adv_blocks;

 for (auto const& bl : gf_struct) {
  block_names.push_back(bl.first);
  int n = bl.second.size();

  index_visitor iv;
  for (auto & ind: bl.second) { apply_visitor(iv, ind); }
  std::vector<std::vector<std::string>> indices{{iv.indices,iv.indices}};

  t_mesh = gf_mesh<retime>{time_window.first, time_window.second, n_t};

  g_ret_blocks.push_back(g_2t_t{{t_mesh, t_mesh}, {n, n}, indices});
  g_adv_blocks.push_back(g_2t_t{{t_mesh, t_mesh}, {n, n}, indices});
 }

 g_ret = make_block_gf(block_names, g_ret_blocks);
 g_adv = make_block_gf(block_names, g_adv_blocks);
}

operator_t try_reduce_to_constant(operator_t const& op, gf_mesh<retime> const& mesh) {
 operator_t res;
 for(auto const& m : op) {
  auto new_coef = try_reduce_to_constant(m.coef, mesh);
  operator_t new_monomial(new_coef);
  for(auto const& c : m.monomial)
   new_monomial *= operator_t::make_canonical(c.stat, c.dagger, c.indices);
  res += new_monomial;
 }
 return res;
}

void solver::solve(solve_parameters_t const& p) {

 // Save parameters
 params = p;

 // Scan for the coefficients constant at all points of t_mesh,
 // and replace them with the corresponding constants
 auto h = try_reduce_to_constant(p.h, t_mesh);

 // Determine fundamental operator set to use
 auto fops = merge(h.make_fundamental_operator_set(), p.h0.make_fundamental_operator_set());

 // Translate bits_per_boson from a dictionary to a list
 using realevol::hilbert_space::Fermion;    // FIXME
 using realevol::hilbert_space::Boson;      // FIXME
 std::vector<int> bits_per_boson(fops.size(Boson), 0);
 for(auto const& b : p.bits_per_boson) bits_per_boson[fops.pos(b.first, Boson)] = b.second;

 auto zero_it = std::find(std::cbegin(bits_per_boson), std::cend(bits_per_boson), 0);
 if(zero_it != bits_per_boson.end())
  TRIQS_RUNTIME_ERROR << "Bosonic index (" << fops.reverse_map(Boson)[*zero_it]
                      << ") has no associated record in bits_per_boson";

 // Full Hilbert space
 class hilbert_space full_hs(fops, bits_per_boson);

 // Analyse structure of the Hilbert space
 hilbert_space_structure hs_struct(h, fops, full_hs, true, t_mesh);


}

}
