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
#include <triqs/test_tools/arrays.hpp>

#include <random>

#include <realevol/sort2.hpp>

using namespace realevol;

TEST(sort2, Random) {

 std::vector<int> v1(100), v2(100);
 // Fill v1 with a random permutation of numbers 0 ... 99
 std::iota(v1.begin(), v1.end(), 0);
 std::mt19937 g(77777);
 std::shuffle(v1.begin(), v1.end(), g);
 ASSERT_FALSE(std::is_sorted(v1.begin(), v1.end()));

 // v2[i] = v1[i] + 10
 std::transform(v1.begin(), v1.end(), v2.begin(), [](int x) { return x + 10; });

 // Sorting
 auto p = sort_permutation(v1);
 auto v3 = apply_permutation(v2, p);
 EXPECT_TRUE(std::is_sorted(v3.begin(), v3.end()));
 apply_permutation_in_place(v2, p);
 EXPECT_TRUE(std::is_sorted(v2.begin(), v2.end()));
}

MAKE_MAIN;
