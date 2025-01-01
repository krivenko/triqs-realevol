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

#include <triqs/utility/first_include.hpp>

// clang-format off
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
// clang-format on

#include <triqs/gfs.hpp>
#include <triqs/mesh/bases/segment.hpp>

#include <realevol/interpolator.hpp>

using namespace realevol;

using triqs::mesh::segment_mesh;

TEST(interpolator, out_of_bounds) {
  segment_mesh m(0,0.5,6);
  interpolator1d<double> interp(m, {3,4,5,6,7,8});

  EXPECT_THROW(interp(-1), triqs::exception);
  EXPECT_THROW(interp(0.6), triqs::exception);
}

TEST(interpolator, evaluate) {
  segment_mesh m(0,0.5,6);
  interpolator1d<double> interp(m, {0,1,4,9,16,25});

  std::vector<double> ref = {0,0.5,1,2.5,4,6.5,9,12.5,16,20.5,25};
  for(int i = 0; i < ref.size(); ++i) {
    EXPECT_CLOSE(ref[i], interp(i*0.05));
  }
}

TEST(interpolator, is_constant) {
  segment_mesh m(0,0.5,6);

  interpolator1d<double> interp1(m, {2,2,2,2,2,2});
  EXPECT_TRUE(is_constant(interp1));
  EXPECT_FALSE(interp1.is_zero());

  interpolator1d<double> interp2(m, 3);
  EXPECT_TRUE(is_constant(interp2));
  EXPECT_FALSE(interp2.is_zero());

  interpolator1d<double> interp3(m, 0);
  EXPECT_TRUE(is_constant(interp3));
  EXPECT_TRUE(interp3.is_zero());
}
