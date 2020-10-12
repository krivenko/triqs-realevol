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
#include <string>
#include <vector>

#include "utility.hpp"

namespace realevol {

std::pair<block_gf_2t_t,block_gf_2t_t>
make_gf_ret_adv(block_gf_2t_t const& g_l, block_gf_2t_t const& g_g) {
 if(g_l.mesh() != g_g.mesh()) TRIQS_RUNTIME_ERROR
  << "g_l and g_g have different blocks.";

 std::vector<gf_2t_t> blocks;
 for(auto bl : g_l.mesh()) {
  auto const& mesh = g_l[bl].mesh();
  if(mesh != g_g[bl].mesh()) TRIQS_RUNTIME_ERROR
   << "Block " << g_l.name[bl]
   << " has different meshes within g_l and g_g.";
  auto shape = get_target_shape(g_l[bl]);
  if(shape != get_target_shape(g_g[bl])) TRIQS_RUNTIME_ERROR
   << "Block " << g_l.name[bl]
   << " has different target shapes within g_l and g_g.";
  auto indices = g_l[bl].indices();

  blocks.push_back(gf_2t_t{mesh, shape, indices, g_l.name});
 }

 auto res = std::make_pair(make_block_gf(g_l.domain().names(), blocks),
                           make_block_gf(g_l.domain().names(), blocks));

 gf_mesh<retime>::mesh_point_t t, tp;
 for(auto bl : g_l.mesh()) {
  auto const& g_l_block = g_l[bl];
  auto const& g_g_block = g_g[bl];
  auto & g_ret_block = res.first[bl];
  auto & g_adv_block = res.second[bl];

  for(auto ttp : g_l_block.mesh()) {
   std::tie(t, tp) = ttp.components_tuple();
   if(t >= tp) g_ret_block[ttp] = g_g_block[ttp] - g_l_block[ttp];
   if(t <= tp) g_adv_block[ttp] = g_l_block[ttp] - g_g_block[ttp];
  }
 }

 return res;
}

}
