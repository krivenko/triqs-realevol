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

#include <set>

#include <triqs/test_tools/arrays.hpp>

#include <realevol/time_expr.hpp>

#include <realevol/hs_structure.hpp>
#include <realevol/mesh_utils.hpp>

using namespace realevol;
using namespace realevol::operators;
using namespace realevol::hilbert_space;

using state_set_t = std::set<fock_state_t>;

TEST(hs_structure, BosonFermion) {

 double U = 3.0;
 double mu = U/2;
 double Omega = 0.5;

 auto h = -mu*(n<time_expr>("up",0) + n<time_expr>("dn",0));
 h += U*n<time_expr>("up",0)*n<time_expr>("dn",0);
 h += Omega*a_dag<time_expr>("B",0)*a<time_expr>("B",0);

 fundamental_operator_set fops;
 fops.insert_fermion("dn",0);
 fops.insert_fermion("up",0);
 fops.insert_boson("B",0);

 triqs::gfs::segment_mesh mesh(0,1.0,101);
 class hilbert_space full_hs(fops, {3});
 hilbert_space_structure<time_expr_operator_t> hss(h, fops, full_hs, fops,
                             is_zero_on_mesh<triqs::gfs::segment_mesh>(mesh));

 EXPECT_EQ(32, hss.sub_hilbert_spaces.size());

 std::set<state_set_t> part, part_ref;
 for(auto const& sp : hss.sub_hilbert_spaces) {
  auto all_states = sp.get_all_fock_states();
  part.insert(state_set_t(all_states.begin(), all_states.end()));
 }
 for(fock_state_t i = 0; i < 32; ++i) part_ref.insert(state_set_t{i});
 EXPECT_EQ(part_ref, part);
}

TEST(hs_structure, BosonFermionCoupled) {

 double U = 3.0;
 double mu = U/2;
 double Omega = 0.5;
 time_expr lambda = "0.2*cos(t)";

 auto h = -mu*(n<time_expr>("up",0) + n<time_expr>("dn",0));
 h += U*n<time_expr>("up",0)*n<time_expr>("dn",0);
 h += Omega*a_dag<time_expr>("B",0)*a<time_expr>("B",0);
 h += lambda*0.5*(n<time_expr>("up",0) - n<time_expr>("dn",0))
            *(a_dag<time_expr>("B",0) + a<time_expr>("B",0));

 fundamental_operator_set fops;
 fops.insert_fermion("dn",0);
 fops.insert_fermion("up",0);
 fops.insert_boson("B",0);

 triqs::gfs::segment_mesh mesh(0,1.0,101);
 class hilbert_space full_hs(fops, {3});
 hilbert_space_structure<time_expr_operator_t> hss(h, fops, full_hs, fops,
                             is_zero_on_mesh<triqs::gfs::segment_mesh>(mesh));

 EXPECT_EQ(4, hss.sub_hilbert_spaces.size());

 std::set<state_set_t> part, part_ref;
 for(auto const& sp : hss.sub_hilbert_spaces) {
  auto all_states = sp.get_all_fock_states();
  part.insert(state_set_t(all_states.begin(), all_states.end()));
 }

 for(fock_state_t f_state : {0,1,2,3}) {
  state_set_t f_sector;
  for(fock_state_t b_state = 0; b_state < 8; ++b_state)
   f_sector.insert(f_state + (b_state << 2));
  part_ref.insert(f_sector);
 }

 EXPECT_EQ(part_ref, part);
}

MAKE_MAIN;
