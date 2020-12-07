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

#include <triqs/gfs.hpp>

#include <triqs/test_tools/arrays.hpp>

#include <realevol/time_interp.hpp>
#include <realevol/mesh_utils.hpp>

using namespace realevol;
using namespace triqs::arrays;
using triqs::gfs::segment_mesh;

segment_mesh m(0, 10, 11);

template<typename F>
auto init_array(segment_mesh const& m_, F && f) -> triqs::arrays::array<decltype(f(.0)), 1> {
  triqs::arrays::array<decltype(f(.0)), 1> res(m_.size());
  for(int i : range(m_.size())) res(i) = f(m_[i]);
  return res;
}

time_interp ti0(m);
time_interp ti1(m, init_array(m, [](double t){ return t*t; }));
time_interp ti2(m, init_array(m, [](double t){ return t + 1; }));
time_interp ti3(m, 4.5);
time_interp ti4(m, init_array(m, [](double t){ return t*t; }), 1.0);
time_interp ti5(m, init_array(m, [](double t){ return t + 1; }),
                   init_array(m, [](double t){ return t*t; })
);
time_interp ti6(m, 4.5, 2.9);
time_interp ti7(m, 1.9+2.8i);
time_interp ti8(m, init_array(m, [](double t){ return t*t + 1i*(t + 1); }));

constexpr int N_t = 5;
double T[N_t] = {0, 0.1, 1, 5.5, 10};
std::complex<double> TI0_res[] = {0, 0, 0, 0, 0};
std::complex<double> TI1_res[] = {0, 0.1, 1, 30.5, 100};
std::complex<double> TI2_res[] = {1, 1.1, 2, 6.5, 11};
std::complex<double> TI3_res[] = {4.5, 4.5, 4.5, 4.5, 4.5};
std::complex<double> TI4_res[] = {{0,1.0}, {0.1,1.0}, {1,1.0}, {30.5,1.0}, {100,1.0}};
std::complex<double> TI5_res[] = {{1,0}, {1.1,0.1}, {2,1}, {6.5,30.5}, {11,100}};
std::complex<double> TI6_res[] = {{4.5,2.9}, {4.5,2.9}, {4.5,2.9}, {4.5,2.9}, {4.5,2.9}};
std::complex<double> TI7_res[] = {{1.9,2.8}, {1.9,2.8}, {1.9,2.8}, {1.9,2.8}, {1.9,2.8}};
std::complex<double> TI8_res[] = {{0,1}, {0.1,1.1}, {1,2}, {30.5,6.5}, {100,11}};

TEST(time_interp, Assignment) {
  time_interp tia(m);
  tia = ti1; EXPECT_EQ(ti1, tia);
}

TEST(time_interp, is_constant) {
  EXPECT_TRUE(is_constant(ti0));
  EXPECT_FALSE(is_constant(ti1));
  EXPECT_FALSE(is_constant(ti2));
  EXPECT_TRUE(is_constant(ti3));
  EXPECT_FALSE(is_constant(ti4));
  EXPECT_FALSE(is_constant(ti5));
  EXPECT_TRUE(is_constant(ti6));
  EXPECT_TRUE(is_constant(ti7));
  EXPECT_FALSE(is_constant(ti8));
}

TEST(time_interp, RealImag) {

  using namespace triqs::utility;
  EXPECT_TRUE(is_zero(ti0));
  EXPECT_FALSE(is_zero(ti1));
  EXPECT_EQ(time_interp(m, 4.5), ti3.real());
  EXPECT_EQ(time_interp(m), ti3.imag());
  EXPECT_EQ(time_interp(m, 4.5), conj(ti3));
  auto ti5_real_ref = time_interp(m, init_array(m, [](double t){ return t + 1; }));
  EXPECT_EQ(ti5_real_ref, ti5.real());
  auto ti5_imag_ref = time_interp(m, init_array(m, [](double t){ return t*t; }));
  EXPECT_EQ(ti5_imag_ref, ti5.imag());
  auto ti5_conj_ref = time_interp(m, init_array(m, [](double t){ return t + 1 - 1i*t*t; }));
  EXPECT_EQ(ti5_conj_ref, conj(ti5));
}

TEST(time_interp, Evaluation) {
  for(int i = 0; i < N_t; ++i) {
    double t = T[i];
    EXPECT_CLOSE(TI0_res[i], ti0(t));
    EXPECT_CLOSE(TI1_res[i], ti1(t));
    EXPECT_CLOSE(TI2_res[i], ti2(t));
    EXPECT_CLOSE(TI3_res[i], ti3(t));
    EXPECT_CLOSE(TI4_res[i], ti4(t));
    EXPECT_CLOSE(TI5_res[i], ti5(t));
    EXPECT_CLOSE(TI6_res[i], ti6(t));
    EXPECT_CLOSE(TI7_res[i], ti7(t));
    EXPECT_CLOSE(TI8_res[i], ti8(t));
  }
}

TEST(time_interp, UnaryMinus) {
  time_interp mti2 = -ti2;
  time_interp mti5 = -ti5;
  for(int i = 0; i < N_t; ++i) {
    double t = T[i];
    EXPECT_CLOSE(-TI2_res[i], mti2(t));
    EXPECT_CLOSE(-TI5_res[i], mti5(t));
  }
}

TEST(time_interp, Addition) {
  time_interp ti1pti2 = ti1 + ti2;
  time_interp ti1phalf = ti1 + 0.5;
  time_interp halfpti2 = 0.5 + ti2;
  time_interp ti1pihalf = ti1 + 0.5i;
  time_interp ihalfpti2 = 0.5i + ti2;

  for(int i = 0; i < N_t; ++i) {
    double t = T[i];
    EXPECT_CLOSE(TI1_res[i]+TI2_res[i], ti1pti2(t));
    EXPECT_CLOSE(TI1_res[i]+0.5, ti1phalf(t));
    EXPECT_CLOSE(0.5+TI2_res[i], halfpti2(t));
    EXPECT_CLOSE(TI1_res[i]+0.5i, ti1pihalf(t));
    EXPECT_CLOSE(0.5i+TI2_res[i], ihalfpti2(t));
  }
}

TEST(time_interp, Subtraction) {
  time_interp ti1mti2 = ti1 - ti2;
  time_interp ti1mhalf = ti1 - 0.5;
  time_interp halfmti2 = 0.5 - ti2;
  time_interp ti1mihalf = ti1 - 0.5i;
  time_interp ihalfmti2 = 0.5i - ti2;

  for(int i = 0; i < N_t; ++i) {
    double t = T[i];
    EXPECT_CLOSE(TI1_res[i]-TI2_res[i], ti1mti2(t));
    EXPECT_CLOSE(TI1_res[i]-0.5, ti1mhalf(t));
    EXPECT_CLOSE(0.5-TI2_res[i], halfmti2(t));
    EXPECT_CLOSE(TI1_res[i]-0.5i, ti1mihalf(t));
    EXPECT_CLOSE(0.5i-TI2_res[i], ihalfmti2(t));
  }
}

TEST(time_interp, Multiplication) {
  time_interp ti1ppti2 = ti1 * ti2;
  time_interp ti5ppti6 = ti5 * ti6;
  time_interp ti1pphalf = ti1 * 0.5;
  time_interp halfppti2 = 0.5 * ti2;
  time_interp ti1ppihalf = ti1 * 0.5i;
  time_interp ihalfppti2 = 0.5i * ti2;

  std::complex<double> TI1_TI2_res[] = {0, 0.2, 2, 201, 1100};
  std::complex<double> TI5_TI6_res[] = {{4.5,2.9}, {4.66,3.64}, {6.1,10.3}, {-59.2,156.1}, {-240.5,481.9}};

  for(long i = 0; i < N_t; ++i) {
    double t = T[i];
    EXPECT_CLOSE(TI1_TI2_res[i], ti1ppti2(t));
    EXPECT_CLOSE(TI5_TI6_res[i], ti5ppti6(t));
    EXPECT_CLOSE(TI1_res[i]*0.5, ti1pphalf(t));
    EXPECT_CLOSE(0.5*TI2_res[i], halfppti2(t));
    EXPECT_CLOSE(TI1_res[i]*0.5i, ti1ppihalf(t));
    EXPECT_CLOSE(0.5i*TI2_res[i], ihalfppti2(t));
  }
}

TEST(time_interp, Division) {
  time_interp ti1dti2 = ti1 / ti2;
  time_interp ti5dti6 = ti5 / ti6;
  time_interp ti1dhalf = ti1 / 0.5;
  time_interp halfdti2 = 0.5 / ti2;
  time_interp ti1dihalf = ti1 / 0.5i;
  time_interp ihalfdti2 = 0.5i / ti2;

  std::complex<double> TI1_TI2_res[] = {0, 0.05, 0.5, 0.5*(25/6.0 + 36/7.0), 100.0 / 11};
  std::complex<double> halfd_TI2_res[] = {0.5, 0.475, 0.25, 0.25*(1/6.0 + 1/7.0), 0.5 / 11};

  for(int i = 0; i < N_t; ++i){
    double t = T[i];
    EXPECT_CLOSE(TI1_TI2_res[i], ti1dti2(t));
    EXPECT_CLOSE(TI5_res[i] / TI6_res[i], ti5dti6(t));
    EXPECT_CLOSE(TI1_res[i] / 0.5, ti1dhalf(t));
    EXPECT_CLOSE(halfd_TI2_res[i], halfdti2(t));
    EXPECT_CLOSE(TI1_res[i] / 0.5i, ti1dihalf(t));
    EXPECT_CLOSE(1i*halfd_TI2_res[i], ihalfdti2(t));
  }
}

TEST(time_interp, try_reduce_to_constant) {

  segment_mesh m2(0,100,101);
  time_interp ti(m2, init_array(m2, [](double t){ return double(t > 50); }));

  EXPECT_FALSE(is_constant(ti));

  segment_mesh m_low(0,50,51), m_mid(25,75,51), m_high(51,100,50);

  auto ti_low = try_reduce_to_constant(ti, m_low);
  auto ti_mid = try_reduce_to_constant(ti, m_mid);
  auto ti_high = try_reduce_to_constant(ti, m_high);

  EXPECT_TRUE(is_constant(ti_low));
  EXPECT_FALSE(is_constant(ti_mid));
  EXPECT_TRUE(is_constant(ti_high));
}

MAKE_MAIN;
