/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2024, I. Krivenko, M. Danilov, P. Kubiczek
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
#include <cmath>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <realevol/time_expr.hpp>
#include <realevol/init_state.hpp>

using namespace realevol;
using namespace realevol::operators;
using namespace realevol::hilbert_space;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  ASSERT(mpi::has_env);
  mpi::environment env(argc, argv);
  return RUN_ALL_TESTS();
}

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

 auto h0 = -mu*(n("up") + n("dn"));
 h0 += eta*(c_dag("up")*c("dn") + c_dag("dn")*c("up"));
 h0 += U*n("up")*n("dn");
 h0 += Omega*a_dag("B")*a("B");
 h0 += lambda*(n("up") + n("dn"))*(a_dag("B") + a("B"));

 eq_solver_parameters_t params;
 params.verbosity = 2;
 params.arpack_min_matrix_size = 257;

 auto ist = make_equilibrium_init_state(h0, fops, T, params, {{{"B"},8}});

 //h5::file f("init_state_equilibrium.ref.h5", 'r');
 //h5_write(f, "real", ist);
 h5::file f("init_state_equilibrium.ref.h5", 'r');
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
  EXPECT_CLOSE(1.0, nda::blas::dot(a, a));
  EXPECT_CLOSE(1.0, abs(nda::blas::dot(wst_ref[i].state.amplitudes(), a)));
 }
 EXPECT_CLOSE(1.0, total_weight);
}

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

 auto h0 = -mu*(n("up") + n("dn"));
 h0 += eta*1i*(c_dag("up")*c("dn") - c_dag("dn")*c("up"));
 h0 += U*n("up")*n("dn");
 h0 += Omega*a_dag("B")*a("B");
 h0 += lambda*1i*(n("up") + n("dn"))*(a_dag("B") - a("B"));

 eq_solver_parameters_t params;
 params.verbosity = 2;
 params.arpack_min_matrix_size = 129;

 auto ist = make_equilibrium_init_state(h0, fops, T, params, {{{"B"},7}});

 //h5::file f("init_state_equilibrium.ref.h5", 'w');
 //h5_write(f, "complex", ist);
 h5::file f("init_state_equilibrium.ref.h5", 'r');
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
  EXPECT_CLOSE(1.0, nda::blas::dotc(a,a));
  EXPECT_CLOSE(1.0, abs(nda::blas::dotc(wst_ref[i].state.amplitudes(),a)));
 }
 EXPECT_CLOSE(1.0, total_weight);
}
