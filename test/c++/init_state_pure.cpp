/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <triqs/test_tools/arrays.hpp>

#include <realevol/time_expr.hpp>
#include <realevol/init_state.hpp>

using namespace realevol;
using namespace realevol::operators;
using namespace realevol::hilbert_space;

TEST(init_state, pure) {
 fundamental_operator_set fops;
 fops.insert_fermion("dn");
 fops.insert_fermion("up");
 fops.insert_boson("B");

 auto generator = (1.0 + c_dag<time_expr>("up") + c_dag<time_expr>("dn")
                       + c_dag<time_expr>("dn")*c_dag<time_expr>("up"))
                      *a_dag<time_expr>("B")*a_dag<time_expr>("B");

 auto st = make_pure_init_state(generator, fops, {{{"B"},3}});

 auto const& wst = st.get_weighted_states()[0];
 auto const& hs = wst.state.get_hilbert();

 EXPECT_EQ(1, st.size());
 EXPECT_EQ(1.0, wst.weight);

 {
  auto f = st.get_full_hs().get_fock_state(fops, {}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.5, wst.state(hs.get_state_index(f)));
 }
 {
  auto f = st.get_full_hs().get_fock_state(fops, {{"dn"}}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.5, wst.state(hs.get_state_index(f)));
 }
 {
  auto f = st.get_full_hs().get_fock_state(fops, {{"up"}}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.5, wst.state(hs.get_state_index(f)));
 }
 {
  auto f = st.get_full_hs().get_fock_state(fops, {{"dn"},{"up"}}, {{"B"},{"B"}});
  EXPECT_CLOSE(0.5, wst.state(hs.get_state_index(f)));
 }

 {
  h5::file ff("init_state_pure.h5", 'w');
  h5_write(ff, "init_state", st);
 }
 {
  h5::file ff("init_state_pure.h5", 'r');
  init_state st_new;
  h5_read(ff, "init_state", st_new);

  EXPECT_EQ(st.size(), st_new.size());
  EXPECT_EQ(st.get_fops(), st_new.get_fops());
  EXPECT_EQ(st.get_full_hs(), st_new.get_full_hs());
  EXPECT_EQ(st.get_sub_hilbert_spaces(), st_new.get_sub_hilbert_spaces());
  auto const& wst_new = st.get_weighted_states()[0];
  EXPECT_EQ(wst.weight, wst_new.weight);
  EXPECT_EQ(wst.state.get_hilbert(), wst_new.state.get_hilbert());
  EXPECT_ARRAY_EQ(wst.state.amplitudes(), wst_new.state.amplitudes());
 }
}

MAKE_MAIN;
