#include <triqs/test_tools/arrays.hpp>

#include <cmath>

#include "init_state.hpp"

using namespace realevol;
//using namespace triqs::operators;
//using namespace triqs::hilbert_space;
using namespace realevol::operators;       // FIXME
using namespace realevol::hilbert_space;   // FIXME

TEST(init_state_equilibrium, real) {
 double U = 3.0;
 double mu = U/2;
 double eta = 0.01;
 double Omega = 0.7;
 double lambda = 0.2;
 double T = 0.1; // beta = 10

 fundamental_operator_set fops;
 fops.insert_fermion("dn");
 fops.insert_fermion("up");
 fops.insert_boson("B");

 auto h0 = -mu*(n<time_expr>("up") + n<time_expr>("dn"));
 h0 += eta*(c_dag<time_expr>("up")*c<time_expr>("dn") + c_dag<time_expr>("dn")*c<time_expr>("up"));
 h0 += U*n<time_expr>("up")*n<time_expr>("dn");
 h0 += Omega*a_dag<time_expr>("B")*a<time_expr>("B");
 h0 += lambda*(n<time_expr>("up") + n<time_expr>("dn"))
             *(a_dag<time_expr>("B") + a<time_expr>("B"));

 eq_solver_parameters_t params;
 params.verbosity = 2;
 params.arpack_min_matrix_size = 257;

 auto ist = make_equilibrium_init_state(h0, fops, T, params, {{{"B"},8}});

 //h5::file f("init_state_equilibrium.ref.h5", H5F_ACC_TRUNC);
 //h5_write(f, "real", ist);
 h5::file f("init_state_equilibrium.ref.h5", H5F_ACC_RDONLY);
 init_state ist_ref;
 h5_read(f, "real", ist_ref);

 EXPECT_EQ(ist_ref.get_fops(), ist.get_fops());
 EXPECT_EQ(ist_ref.get_full_hs(), ist.get_full_hs());

 auto const& hs = ist.get_sub_hilbert_spaces();
 auto const& hs_ref = ist_ref.get_sub_hilbert_spaces();
 ASSERT_EQ(hs_ref.size(), hs.size());
 for(int i = 0; i < hs.size(); ++i) EXPECT_EQ(hs_ref[i],hs[i]);

 auto const& wst = ist.get_weighted_states();
 auto const& wst_ref = ist_ref.get_weighted_states();
 double total_weight = 0;
 ASSERT_EQ(wst_ref.size(), wst.size());
 for(int i = 0; i < wst.size(); ++i) {
  EXPECT_CLOSE(wst_ref[i].weight, wst[i].weight);
  total_weight += wst[i].weight;
  EXPECT_EQ(wst_ref[i].state.get_hilbert().get_index(),
            wst[i].state.get_hilbert().get_index());
  auto const& a = wst[i].state.amplitudes();
  EXPECT_CLOSE(1.0, dot(a,a));
  EXPECT_CLOSE(1.0, abs(dot(wst_ref[i].state.amplitudes(),a)));
 }
 EXPECT_CLOSE(1.0, total_weight);
}
/*
TEST(init_state_equilibrium, complex) {
 double U = 3.0;
 double mu = U/2;
 double eta = 0.01;
 double Omega = 0.7;
 double lambda = 0.2;
 double T = 0.1;

 fundamental_operator_set fops;
 fops.insert_fermion("dn");
 fops.insert_fermion("up");
 fops.insert_boson("B");

 auto h0 = -mu*(n<time_expr>("up") + n<time_expr>("dn"));
 h0 += eta*1_j*(c_dag<time_expr>("up")*c<time_expr>("dn") - c_dag<time_expr>("dn")*c<time_expr>("up"));
 h0 += U*n<time_expr>("up")*n<time_expr>("dn");
 h0 += Omega*a_dag<time_expr>("B")*a<time_expr>("B");
 h0 += lambda*1_j*(n<time_expr>("up") + n<time_expr>("dn"))
             *(a_dag<time_expr>("B") - a<time_expr>("B"));

 eq_solver_parameters_t params;
 params.verbosity = 2;
 params.arpack_min_matrix_size = 129;

 auto ist = make_equilibrium_init_state(h0, fops, T, params, {{{"B"},7}});

 //h5::file f("init_state_equilibrium.ref.h5", H5F_ACC_RDWR);
 //h5_write(f, "complex", ist);
 h5::file f("init_state_equilibrium.ref.h5", H5F_ACC_RDONLY);
 init_state ist_ref;
 h5_read(f, "complex", ist_ref);

 EXPECT_EQ(ist_ref.get_fops(), ist.get_fops());
 EXPECT_EQ(ist_ref.get_full_hs(), ist.get_full_hs());

 auto const& hs = ist.get_sub_hilbert_spaces();
 auto const& hs_ref = ist_ref.get_sub_hilbert_spaces();
 ASSERT_EQ(hs_ref.size(), hs.size());
 for(int i = 0; i < hs.size(); ++i) EXPECT_EQ(hs_ref[i],hs[i]);

 auto const& wst = ist.get_weighted_states();
 auto const& wst_ref = ist_ref.get_weighted_states();
 double total_weight = 0;
 ASSERT_EQ(wst_ref.size(), wst.size());
 for(int i = 0; i < wst.size(); ++i) {
  EXPECT_CLOSE(wst_ref[i].weight, wst[i].weight);
  total_weight += wst[i].weight;
  EXPECT_EQ(wst_ref[i].state.get_hilbert().get_index(),
            wst[i].state.get_hilbert().get_index());
  auto const& a = wst[i].state.amplitudes();
  EXPECT_CLOSE(1.0, dotc(a,a));
  EXPECT_CLOSE(1.0, abs(dotc(wst_ref[i].state.amplitudes(),a)));
 }
 EXPECT_CLOSE(1.0, total_weight);
}
*/
MAKE_MAIN;
