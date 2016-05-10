#include <triqs/test_tools/arrays.hpp>

#include <cmath>

#include "init_state.hpp"

using namespace realevol;
using operators::c;
using operators::c_dag;
using operators::n;
using operators::a;
using operators::a_dag;
using namespace hilbert_space;

TEST(init_state, pure) {

 double U = 3.0;
 double mu = U/2;
 double Omega = 0.5;
 time_expr lambda = "0.2*cos(t)";

 auto h = -mu*(n<time_expr>("up") + n<time_expr>("dn"));
 h += U*n<time_expr>("up")*n<time_expr>("dn");
 h += Omega*a_dag<time_expr>("B")*a<time_expr>("B");
 h += lambda*0.5*(n<time_expr>("up") - n<time_expr>("dn"))
            *(a_dag<time_expr>("B") + a<time_expr>("B"));

 fundamental_operator_set fops;
 fops.insert_fermion("dn");
 fops.insert_fermion("up");
 fops.insert_boson("B");

 triqs::gfs::segment_mesh mesh(0,1.0,101);
 hilbert_space_structure hss(fops, h, {3}, true, mesh);

 auto generator = 0.25*(1.0 + c_dag<time_expr>("up") + c_dag<time_expr>("dn")
                            + c_dag<time_expr>("dn")*c_dag<time_expr>("up"))
                      *a_dag<time_expr>("B")*a_dag<time_expr>("B");

 auto st = init_state_pure(generator, hss);

 EXPECT_EQ(1, st.size());
 EXPECT_EQ(1.0, st[0].weight);
 EXPECT_EQ(4, st[0].parts.size());

 {
  auto f = hss.full_hs.get_fock_state(fops, {}, {{"B"},{"B"}});
  auto hs = st[0].parts[0].get_hilbert();
  EXPECT_CLOSE(0.25*std::sqrt(2.0), st[0].parts[0](hs.get_state_index(f)));
 }
 {
  auto f = hss.full_hs.get_fock_state(fops, {{"dn"}}, {{"B"},{"B"}});
  auto hs = st[0].parts[1].get_hilbert();
  EXPECT_CLOSE(0.25*std::sqrt(2.0), st[0].parts[1](hs.get_state_index(f)));
 }
 {
  auto f = hss.full_hs.get_fock_state(fops, {{"up"}}, {{"B"},{"B"}});
  auto hs = st[0].parts[2].get_hilbert();
  EXPECT_CLOSE(0.25*std::sqrt(2.0), st[0].parts[2](hs.get_state_index(f)));
 }
 {
  auto f = hss.full_hs.get_fock_state(fops, {{"dn"},{"up"}}, {{"B"},{"B"}});
  auto hs = st[0].parts[3].get_hilbert();
  EXPECT_CLOSE(0.25*std::sqrt(2.0), st[0].parts[3](hs.get_state_index(f)));
 }
}

TEST(init_state, thermal) {
}

MAKE_MAIN;
