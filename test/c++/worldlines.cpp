/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2023, I. Krivenko, M. Danilov, P. Kubiczek
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
#include <numeric>
#include <vector>

#include <mpi/mpi.hpp>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <triqs/mesh/retime.hpp>

#include <realevol/array_utility.hpp>
#include <realevol/time_expr.hpp>
#include <realevol/make_static_op.hpp>
#include <realevol/mesh_utils.hpp>
#include <realevol/init_state.hpp>
#include <realevol/worldlines.hpp>

using namespace realevol;
using namespace realevol::operators;
using namespace realevol::hilbert_space;

namespace realevol {

using time_expr_operator_t = realevol::operators::many_body_operator_generic<time_expr>;

template<std::size_t NPoints>
bool operator==(worldline_desc_t<NPoints> const& wl1,
                worldline_desc_t<NPoints> const& wl2) {
  return wl1.factor == wl2.factor &&
          wl1.M == wl2.M &&
          wl1.sp_indices == wl2.sp_indices &&
          wl1.weighted_state_index == wl2.weighted_state_index;
}

} // namespace realevol

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  ASSERT(mpi::has_env);
  mpi::environment env(argc, argv);
  return RUN_ALL_TESTS();
}

class WorldlinesTest : public ::testing::Test {

protected:

  init_state initial_state;
  std::unique_ptr<hilbert_space_structure<time_expr_operator_t>> hss1, hss2;
  branchings_t branchings1, branchings2;
  std::unique_ptr<worldlines_maker<time_expr_operator_t>> wlm1, wlm2;

  // Model: Single Hubbard atom + spin flips switched on at t = 0
  const double U = 3.0;
  const double mu = 1.6;
  const double h = 0.2;
  const double T = 0.5;
  const double J = 0.1;
  mesh::retime t_mesh = {0, 1.0, 10};

  void SetUp() override {
    fundamental_operator_set fops;
    fops.insert_fermion("dn");
    fops.insert_fermion("up");

    auto H0 = (-mu - h) * n<time_expr>("dn") +
              (-mu + h) * n<time_expr>("up") +
              U * n<time_expr>("up") * n<time_expr>("dn");

    eq_solver_parameters_t params;
    params.verbosity = 0;
    params.arpack_min_matrix_size = 8;

    initial_state = make_equilibrium_init_state(make_static_op(H0), fops, T, params);

    // H1 mixes |up> and |dn>

    time_expr_operator_t H1 = J*(c_dag<time_expr>("dn") * c<time_expr>("up") +
                                 c_dag<time_expr>("up") * c<time_expr>("dn"));

    hss1 = std::make_unique<hilbert_space_structure<time_expr_operator_t>>(
      H1,
      initial_state.get_fops(),
      initial_state.get_full_hs(),
      initial_state.get_fops(),
      is_zero_on_mesh<mesh::retime>(t_mesh)
    );

    branchings1 = hss1->compute_branchings(
      initial_state.get_sub_hilbert_spaces()
    );

    wlm1 = std::make_unique<worldlines_maker<time_expr_operator_t>>(
      initial_state, *hss1, branchings1
    );

    // H2 does not mix any atomic states
    time_expr_operator_t H2;

    hss2 = std::make_unique<hilbert_space_structure<time_expr_operator_t>>(
      H2,
      initial_state.get_fops(),
      initial_state.get_full_hs(),
      initial_state.get_fops(),
      is_zero_on_mesh<mesh::retime>(t_mesh)
    );

    branchings2 = hss2->compute_branchings(
      initial_state.get_sub_hilbert_spaces()
    );

    wlm2 = std::make_unique<worldlines_maker<time_expr_operator_t>>(
      initial_state, *hss2, branchings2
    );
  }

};

TEST_F(WorldlinesTest, InitialState) {
  // Inverse temperature
  double beta = 1/T;

  // Weights of atomic states
  std::vector<double> atomic_weights = {
    1.0,
    std::exp(-beta * (-mu - h)),
    std::exp(-beta * (-mu + h)),
    std::exp(-beta * (-2*mu + U))
  };

  // Atomic partition function
  double Z = std::accumulate(atomic_weights.begin(), atomic_weights.end(), .0);

  // Normalize atomic weights
  std::transform(atomic_weights.begin(),
                 atomic_weights.end(),
                 atomic_weights.begin(),
                 [Z](double w) { return w / Z; }
                );

  auto const& wst = initial_state.get_weighted_states();

  ASSERT_EQ(atomic_weights.size(), wst.size());
  for(int n = 0; n < atomic_weights.size(); ++n) {
    EXPECT_NEAR(atomic_weights[n], wst[n].weight, 1e-10);
    EXPECT_ARRAY_NEAR(nda::vector<dcomplex>{1.0},
                      wst[n].state.amplitudes());
  }
}

TEST_F(WorldlinesTest, Branchings) {
  ASSERT_EQ(hss1->sub_hilbert_spaces.size(), 3);
  // Subspace '1' generated by H1 is spanned by |up> and |dn>
  std::vector<std::set<long>> ref_branchings1 = {{0}, {1}, {1}, {2}};
  ASSERT_EQ(ref_branchings1, branchings1);

  ASSERT_EQ(hss2->sub_hilbert_spaces.size(), 4);
  std::vector<std::set<long>> ref_branchings2 = {{0}, {1}, {2}, {3}};
  ASSERT_EQ(ref_branchings2, branchings2);
}

template<std::size_t NPoints>
worldline_desc_t<NPoints> make_ref_worldline(
  std::complex<double> factor,
  std::array<static_operator_t, NPoints> const& ops,
  std::array<long, NPoints + 1> sp_indices,
  int weighted_state_index
) {
  auto M = map_array<monomial_t>([](auto const& op) {
    return op.begin()->monomial;
  }, ops);
  return worldline_desc_t<NPoints>{
    factor,
    std::move(M),
    std::move(sp_indices),
    weighted_state_index
  };
}

TEST_F(WorldlinesTest, make_worldlines1) {
  auto A = n("up") + 2i * c_dag("up") * c("dn");

  auto wls1 = wlm1->make_worldlines<1>({A});

  std::vector<worldline_desc_t<1>> ref_wls1 = {
    make_ref_worldline<1>(1.0, {n("up")}, {1, 1}, 1),
    make_ref_worldline<1>(1.0, {n("up")}, {1, 1}, 2),
    make_ref_worldline<1>(1.0, {n("up")}, {2, 2}, 3),
    make_ref_worldline<1>(2.0i, {c_dag("up") * c("dn")}, {1, 1}, 1),
    make_ref_worldline<1>(2.0i, {c_dag("up") * c("dn")}, {1, 1}, 2),
    make_ref_worldline<1>(2.0i, {c_dag("up") * c("dn")}, {2, 2}, 3)
  };

  EXPECT_EQ(ref_wls1, wls1);

  auto wls2 = wlm2->make_worldlines<1>({A});

  std::vector<worldline_desc_t<1>> ref_wls2 = {
    make_ref_worldline<1>(1.0, {n("up")}, {2, 2}, 2),
    make_ref_worldline<1>(1.0, {n("up")}, {3, 3}, 3)
  };

  EXPECT_EQ(ref_wls2, wls2);
}


TEST_F(WorldlinesTest, make_worldlines2) {
  auto Sp = c_dag("up") * c("dn");
  auto Sm = c_dag("dn") * c("up");
  auto A = n("up") + 2i * Sp;
  auto B = 3 * n("dn") + 4i * Sm;

  auto wls1 = wlm1->make_worldlines<2>({B, A});

  std::vector<worldline_desc_t<2>> ref_wls1 = {
    make_ref_worldline<2>(4.0i, {Sm, n("up")}, {1, 1, 1}, 1),
    make_ref_worldline<2>(4.0i, {Sm, n("up")}, {1, 1, 1}, 2),
    make_ref_worldline<2>(4.0i, {Sm, n("up")}, {2, 2, 2}, 3),
    make_ref_worldline<2>(-8.0, {Sm, Sp}, {1, 1, 1}, 1),
    make_ref_worldline<2>(-8.0, {Sm, Sp}, {1, 1, 1}, 2),
    make_ref_worldline<2>(-8.0, {Sm, Sp}, {2, 2, 2}, 3),
    make_ref_worldline<2>(3.0, {n("dn"), n("up")}, {1, 1, 1}, 1),
    make_ref_worldline<2>(3.0, {n("dn"), n("up")}, {1, 1, 1}, 2),
    make_ref_worldline<2>(3.0, {n("dn"), n("up")}, {2, 2, 2}, 3),
    make_ref_worldline<2>(6.0i, {n("dn"), Sp}, {1, 1, 1}, 1),
    make_ref_worldline<2>(6.0i, {n("dn"), Sp}, {1, 1, 1}, 2),
    make_ref_worldline<2>(6.0i, {n("dn"), Sp}, {2, 2, 2}, 3)
  };

  EXPECT_EQ(ref_wls1, wls1);

  auto wls2 = wlm2->make_worldlines<2>({B, A});

  std::vector<worldline_desc_t<2>> ref_wls2 = {
    make_ref_worldline<2>(-8.0, {Sm, Sp}, {2, 1, 2}, 2),
    make_ref_worldline<2>(3.0, {n("dn"), n("up")}, {3, 3, 3}, 3),
  };

  EXPECT_EQ(ref_wls2, wls2);
}

TEST_F(WorldlinesTest, make_worldlines3) {
  auto Sp = c_dag("up") * c("dn");
  auto Sm = c_dag("dn") * c("up");
  auto A = 2 * Sp;
  auto B = n("dn") - n("up");
  auto C = 3 * Sm;

  auto wls1 = wlm1->make_worldlines<3>({C, B, A});

  std::vector<worldline_desc_t<3>> ref_wls1 = {
    make_ref_worldline<3>(6.0, {Sm, n("dn"), Sp}, {1, 1, 1, 1}, 1),
    make_ref_worldline<3>(6.0, {Sm, n("dn"), Sp}, {1, 1, 1, 1}, 2),
    make_ref_worldline<3>(6.0, {Sm, n("dn"), Sp}, {2, 2, 2, 2}, 3),
    make_ref_worldline<3>(-6.0, {Sm, n("up"), Sp}, {1, 1, 1, 1}, 1),
    make_ref_worldline<3>(-6.0, {Sm, n("up"), Sp}, {1, 1, 1, 1}, 2),
    make_ref_worldline<3>(-6.0, {Sm, n("up"), Sp}, {2, 2, 2, 2}, 3)
  };

  EXPECT_EQ(ref_wls1, wls1);

  auto wls2 = wlm2->make_worldlines<3>({C, B, A});

  std::vector<worldline_desc_t<3>> ref_wls2 = {
    make_ref_worldline<3>(6.0, {Sm, n("dn"), Sp}, {2, 1, 1, 2}, 2)
  };

  EXPECT_EQ(ref_wls2, wls2);
}
