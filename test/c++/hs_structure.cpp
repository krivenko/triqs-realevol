#include <triqs/test_tools/arrays.hpp>

#include <set>

#include "hs_structure.hpp"

using namespace realevol;
using operators::c;
using operators::c_dag;
using operators::n;
using operators::a;
using operators::a_dag;
using namespace hilbert_space;

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


 hilbert_space_structure hss(fops, h, {3}, true);

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
 hilbert_space_structure hss(fops, h, {3}, true, mesh);

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
