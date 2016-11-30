#include <triqs/test_tools/arrays.hpp>

#include <cmath>

#include "init_state.hpp"

using namespace realevol;
//using namespace triqs::operators;
//using namespace triqs::hilbert_space;
using namespace realevol::operators;       // FIXME
using namespace realevol::hilbert_space;   // FIXME

TEST(init_state, pure) {
 fundamental_operator_set fops;
 fops.insert_fermion("dn");
 fops.insert_fermion("up");
 fops.insert_boson("B");

 auto generator = 0.25*(1.0 + c_dag<time_expr>("up") + c_dag<time_expr>("dn")
                            + c_dag<time_expr>("dn")*c_dag<time_expr>("up"))
                      *a_dag<time_expr>("B")*a_dag<time_expr>("B");

 auto st = make_pure_init_state(generator, fops, {{{"B"},3}});

 auto wst = st.get_weighted_states()[0];
 auto hs = wst.parts[0].get_hilbert();

 EXPECT_EQ(1, st.size());
 EXPECT_EQ(1.0, wst.weight);
 EXPECT_EQ(1, wst.parts.size());

 {
  auto f = st.get_full_hs().get_fock_state(fops, {}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.25*std::sqrt(2.0), wst.parts[0](hs.get_state_index(f)));
 }
 {
  auto f = st.get_full_hs().get_fock_state(fops, {{"dn"}}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.25*std::sqrt(2.0), wst.parts[0](hs.get_state_index(f)));
 }
 {
  auto f = st.get_full_hs().get_fock_state(fops, {{"up"}}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.25*std::sqrt(2.0), wst.parts[0](hs.get_state_index(f)));
 }
 {
  auto f = st.get_full_hs().get_fock_state(fops, {{"dn"},{"up"}}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.25*std::sqrt(2.0), wst.parts[0](hs.get_state_index(f)));
 }
}

MAKE_MAIN;
