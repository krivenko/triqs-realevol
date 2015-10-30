/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2015 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <test_tools.hpp>

#include <cstdlib>
#include <cmath>
#include <complex>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <triqs/gfs/gf.hpp>
#include <triqs/gfs/meshes/segment.hpp>
#include <triqs/utility/complex_ops.hpp>

#include "time_expr.hpp"

using namespace realevol;
using triqs::gfs::segment_mesh;

time_expr te1("t^2"), te2("t + sin(pi/2)"), te3("sqrt(9.0) + 1.5");
time_expr te4("t^2",1.0), te5("t + sin(pi/2)","t^3");
time_expr te6("sqrt(9.0) + 1.5","2.9"), te7(1.9+2.8_j);

double T[] = {0, 0.1, 10, 55};
dcomplex TE1_res[] = {0, 0.01, 100, 3025};
dcomplex TE2_res[] = {1, 1.1, 11, 56};
dcomplex TE3_res[] = {4.5, 4.5, 4.5, 4.5};
dcomplex TE4_res[] = {{0,1.0}, {0.01,1.0}, {100,1.0}, {3025,1.0}};
dcomplex TE5_res[] = {{1,0}, {1.1,0.001}, {11,1000}, {56,166375}};
dcomplex TE6_res[] = {{4.5,2.9}, {4.5,2.9}, {4.5,2.9}, {4.5,2.9}};
dcomplex TE7_res[] = {{1.9,2.8}, {1.9,2.8}, {1.9,2.8}, {1.9,2.8}};

TEST(time_expr,Assignments) {
  time_expr tea;
  tea = te1; EXPECT_EQ(te1,tea);
  tea = "t^2"; EXPECT_EQ(te1,tea);
  tea = std::string("t^2"); EXPECT_EQ(te1,tea);
  tea = 4.5; EXPECT_EQ(te3,tea);
  tea = 0.7_te; EXPECT_EQ(time_expr("0.7"),tea);
  tea = "2*t+4*t^3"_te; EXPECT_EQ(time_expr("2*t+4*t^3"),tea);
  tea = "8-t"_te; EXPECT_EQ(time_expr("8-t"),tea);
  tea = time_expr("t^2",1.0); EXPECT_EQ(time_expr("t^2",1.0),tea);
  tea = time_expr{2.0,5.0}; EXPECT_EQ(time_expr("2.0","5.0"),tea);
  EXPECT_EQ(time_expr(2.0,5.0),tea);
}

TEST(time_expr,Literals) {
  time_expr tel1 = "t^3-2*t";
  time_expr tel2 = 0.3;
  EXPECT_EQ(tel1,"t^3-2*t"_te);
  EXPECT_EQ(tel2,0.3_te);
}

TEST(time_expr,is_constant) {
  EXPECT_FALSE(is_constant(te1));
  EXPECT_FALSE(is_constant(te2));
  EXPECT_TRUE(is_constant(te3));
  EXPECT_FALSE(is_constant(te4));
  EXPECT_FALSE(is_constant(te5));
  EXPECT_TRUE(is_constant(te6));
  EXPECT_TRUE(is_constant(te7));
}

TEST(time_expr,RealImag) {
  time_expr te0 = "0";

  using namespace triqs::utility;
  EXPECT_TRUE(is_zero(te0));
  EXPECT_FALSE(is_zero(te1));
  EXPECT_EQ(time_expr("sqrt(9.0) + 1.5"),te3.real());
  EXPECT_EQ(time_expr(),te3.imag());
  EXPECT_EQ(time_expr("sqrt(9.0) + 1.5"),conj(te3));
  EXPECT_EQ(time_expr("t + sin(pi/2)"),te5.real());
  EXPECT_EQ(time_expr("t^3"),te5.imag());
  EXPECT_EQ(time_expr("t + sin(pi/2)","-(t^3)"),conj(te5));
}

TEST(time_expr,Evaluation) {
  for(int i = 0; i < 4; ++i){
    double t = T[i];
    EXPECT_CLOSE(TE1_res[i],te1(t));
    EXPECT_CLOSE(TE2_res[i],te2(t));
    EXPECT_CLOSE(TE3_res[i],te3(t));
    EXPECT_CLOSE(TE4_res[i],te4(t));
    EXPECT_CLOSE(TE5_res[i],te5(t));
    EXPECT_CLOSE(TE6_res[i],te6(t));
    EXPECT_CLOSE(TE7_res[i],te7(t));
  }
}

TEST(time_expr,UnaryMinus) {
  time_expr mte2 = -te2;
  time_expr mte5 = -te5;
  for(int i = 0; i < 4; ++i) {
    double t = T[i];
    EXPECT_CLOSE(-TE2_res[i],mte2(t));
    EXPECT_CLOSE(-TE5_res[i],mte5(t));
  }
}

TEST(time_expr,Addition) {
  time_expr te1pte2 = te1 + te2;
  time_expr te1phalf = te1 + 0.5;
  time_expr halfpte2 = 0.5 + te2;
  time_expr te1pihalf = te1 + 0.5_j;
  time_expr ihalfpte2 = 0.5_j + te2;

  for(int i = 0; i < 4; ++i){
    double t = T[i];
    EXPECT_CLOSE(TE1_res[i]+TE2_res[i],te1pte2(t));
    EXPECT_CLOSE(TE1_res[i]+0.5,te1phalf(t));
    EXPECT_CLOSE(0.5+TE2_res[i],halfpte2(t));
    EXPECT_CLOSE(TE1_res[i]+0.5_j,te1pihalf(t));
    EXPECT_CLOSE(0.5_j+TE2_res[i],ihalfpte2(t));
  }
}

TEST(time_expr,Subtraction) {
  time_expr te1mte2 = te1 - te2;
  time_expr te1mhalf = te1 - 0.5;
  time_expr halfmte2 = 0.5 - te2;
  time_expr te1mihalf = te1 - 0.5_j;
  time_expr ihalfmte2 = 0.5_j - te2;

  for(int i = 0; i < 4; ++i){
    double t = T[i];
    EXPECT_CLOSE(TE1_res[i]-TE2_res[i],te1mte2(t));
    EXPECT_CLOSE(TE1_res[i]-0.5,te1mhalf(t));
    EXPECT_CLOSE(0.5-TE2_res[i],halfmte2(t));
    EXPECT_CLOSE(TE1_res[i]-0.5_j,te1mihalf(t));
    EXPECT_CLOSE(0.5_j-TE2_res[i],ihalfmte2(t));
  }
}

TEST(time_expr,Multiplication) {
  time_expr te1ppte2 = te1 * te2;
  time_expr te5ppte6 = te5 * te6;
  time_expr te1pphalf = te1 * 0.5;
  time_expr halfppte2 = 0.5 * te2;
  time_expr te1ppihalf = te1 * 0.5_j;
  time_expr ihalfppte2 = 0.5_j * te2;

  for(long i = 0; i < 4; ++i){
    double t = T[i];
    EXPECT_CLOSE(TE1_res[i]*TE2_res[i],te1ppte2(t));
    EXPECT_CLOSE(TE5_res[i]*TE6_res[i],te5ppte6(t));
    EXPECT_CLOSE(TE1_res[i]*0.5,te1pphalf(t));
    EXPECT_CLOSE(0.5*TE2_res[i],halfppte2(t));
    EXPECT_CLOSE(TE1_res[i]*0.5_j,te1ppihalf(t));
    EXPECT_CLOSE(0.5_j*TE2_res[i],ihalfppte2(t));
  }
}

TEST(time_expr,Division) {
  time_expr te1dte2 = te1 / te2;
  time_expr te5dte6 = te5 / te6;
  time_expr te1dhalf = te1 / 0.5;
  time_expr halfdte2 = 0.5 / te2;
  time_expr te1dihalf = te1 / 0.5_j;
  time_expr ihalfdte2 = 0.5_j / te2;

  for(int i = 0; i < 4; ++i){
    double t = T[i];
    EXPECT_CLOSE(TE1_res[i]/TE2_res[i],te1dte2(t));
    EXPECT_CLOSE(TE5_res[i]/TE6_res[i],te5dte6(t));
    EXPECT_CLOSE(TE1_res[i]/0.5,te1dhalf(t));
    EXPECT_CLOSE(0.5/TE2_res[i],halfdte2(t));
    EXPECT_CLOSE(TE1_res[i]/0.5_j,te1dihalf(t));
    EXPECT_CLOSE(0.5_j/TE2_res[i],ihalfdte2(t));
  }
}

TEST(time_expr,try_reduce_to_constant) {
  segment_mesh m(0,100,101);

  time_expr te0t("0*t"), te1t("1*t"),
            te0t2("0*t","2"), te1tt3("1*t","t^3");
  EXPECT_FALSE(is_constant(te0t));
  EXPECT_FALSE(is_constant(te1t));
  EXPECT_FALSE(is_constant(te0t2));
  EXPECT_FALSE(is_constant(te1tt3));

  try_reduce_to_constant(te0t, m);
  try_reduce_to_constant(te1t, m);
  try_reduce_to_constant(te0t2, m);
  try_reduce_to_constant(te1tt3, m);

  EXPECT_TRUE(is_constant(te0t));
  EXPECT_FALSE(is_constant(te1t));
  EXPECT_TRUE(is_constant(te0t2));
  EXPECT_FALSE(is_constant(te1tt3));
}

TEST(time_expr,Serialization) {
  // Write to an archive
  std::stringstream archive_str;
  boost::archive::text_oarchive oa(archive_str);
  oa << te1 << te2 << te3 << te4 << te5 << te6 << te7;

  // Read from an archive
  boost::archive::text_iarchive ia(archive_str);
  time_expr read_expr;
  ia >> read_expr; EXPECT_EQ(te1,read_expr);
  ia >> read_expr; EXPECT_EQ(te2,read_expr);
  ia >> read_expr; EXPECT_EQ(te3,read_expr);
  ia >> read_expr; EXPECT_EQ(te4,read_expr);
  ia >> read_expr; EXPECT_EQ(te5,read_expr);
  ia >> read_expr; EXPECT_EQ(te6,read_expr);
  ia >> read_expr; EXPECT_EQ(te7,read_expr);
}

MAKE_MAIN;
