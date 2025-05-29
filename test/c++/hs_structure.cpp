/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2025, I. Krivenko, M. Danilov, P. Kubiczek
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
#include <set>
#include <vector>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <realevol/time_expr.hpp>

#include <realevol/hs_structure.hpp>
#include <realevol/mesh_utils.hpp>

using namespace realevol;
using namespace realevol::operators;
using namespace realevol::hilbert_space;

using state_set_t = std::set<fock_state_t>;

using time_expr_operator_t = realevol::operators::many_body_operator_generic<time_expr>;

void test_make_monomial_connections(hilbert_space_structure<time_expr_operator_t> const& hss) {

  using realevol::hilbert_space::statistic_enum::Fermion;
  using realevol::hilbert_space::statistic_enum::Boson;
  canonical_ops_t ops[] = {
    canonical_ops_t{Fermion, true, {"dn", 0}},
    canonical_ops_t{Fermion, false, {"dn", 0}},
    canonical_ops_t{Fermion, true, {"up", 0}},
    canonical_ops_t{Fermion, false, {"up", 0}},
    canonical_ops_t{Boson, true, {"B", 0}},
    canonical_ops_t{Boson, false, {"B", 0}}
  };

  auto get_ref_conn = [&hss](bool dagger,
                             realevol::hilbert_space::statistic_enum stat,
                             int linear_index) {
    auto c = (dagger ? hss.creation_connection : hss.annihilation_connection)
             [stat](linear_index, range::all);
    return std::vector<int>(c.begin(), c.end());
  };

  std::vector<int> ref_conns[] = {get_ref_conn(true, Fermion, 0),
                                  get_ref_conn(false, Fermion, 0),
                                  get_ref_conn(true, Fermion, 1),
                                  get_ref_conn(false, Fermion, 1),
                                  get_ref_conn(true, Boson, 0),
                                  get_ref_conn(false, Boson, 0)};

  for(int o = 0; o < 6; ++o) {
    EXPECT_EQ(hss.make_monomial_connections(monomial_t{ops[o]}), ref_conns[o]);
  }

  auto compose = [](std::vector<int> const& v1, std::vector<int> const& v2) {
    std::vector<int> v(v1);
    for(int n = 0; n < v.size(); ++n) {
      if(v[n] == -1) continue;
      v[n] = v2[v[n]];
    }
    return v;
  };

  for(int o1 = 0; o1 < 6; ++o1) {
    for(int o2 = 0; o2 < 6; ++o2) {
      auto conn = hss.make_monomial_connections(monomial_t{ops[o1], ops[o2]});
      EXPECT_EQ(conn, compose(ref_conns[o2], ref_conns[o1]));
    }
  }
}

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

 triqs::mesh::retime mesh(0,1.0,101);
 class hilbert_space full_hs(fops, {3});
 hilbert_space_structure<time_expr_operator_t> hss(h, fops, full_hs, fops,
                             is_zero_on_mesh<triqs::mesh::retime>(mesh));

 EXPECT_EQ(32, hss.sub_hilbert_spaces.size());

 std::set<state_set_t> part, part_ref;
 for(auto const& sp : hss.sub_hilbert_spaces) {
  auto all_states = sp.get_all_fock_states();
  part.insert(state_set_t(all_states.begin(), all_states.end()));
 }
 for(fock_state_t i = 0; i < 32; ++i) part_ref.insert(state_set_t{i});
 EXPECT_EQ(part_ref, part);

 test_make_monomial_connections(hss);
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

 triqs::mesh::retime mesh(0,1.0,101);
 class hilbert_space full_hs(fops, {3});
 hilbert_space_structure<time_expr_operator_t> hss(h, fops, full_hs, fops,
                             is_zero_on_mesh<triqs::mesh::retime>(mesh));

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

 test_make_monomial_connections(hss);
}
