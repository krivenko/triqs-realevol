#include <triqs/test_tools/arrays.hpp>

#include <cmath>

#include "init_state.hpp"

using namespace realevol;
//using namespace triqs::operators;
//using namespace triqs::hilbert_space;
using namespace realevol::operators;       // FIXME
using namespace realevol::hilbert_space;   // FIXME

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
 class hilbert_space full_hs(fops, {3});
 hilbert_space_structure hss(h, fops, full_hs, true, mesh);

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

 double beta = 10;
 double U = 3.0;
 double mu = U/2;
 double Omega = 0.5;
 double lambda = 0.2;

 auto h = -mu*(n<time_expr>("up") + n<time_expr>("dn"));
 h += U*n<time_expr>("up")*n<time_expr>("dn");

 fundamental_operator_set fops;
 fops.insert_fermion("dn");
 fops.insert_fermion("up");
 fops.insert_boson("B");

 triqs::gfs::segment_mesh mesh(0,1.0,101);
 class hilbert_space full_hs(fops, {5});
 hilbert_space_structure hss(h, fops, full_hs, true, mesh);

 auto h0 = -mu*(n<time_expr>("up") + n<time_expr>("dn"));
 h0 += U*n<time_expr>("up")*n<time_expr>("dn");
 h0 += Omega*a_dag<time_expr>("B")*a<time_expr>("B");
 h0 += lambda*(n<time_expr>("up") + n<time_expr>("dn"))
             *(a_dag<time_expr>("B") + a<time_expr>("B"));

 auto st = init_state_thermal(h0, hss, beta);
 // TODO

}

MAKE_MAIN;
