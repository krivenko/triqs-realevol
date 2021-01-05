/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2021, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <map>
#include <vector>

// FIXME: Code in <triqs/utility/variant_extensions.hpp> depends on these
// headers but does not include them.
//
// https://github.com/TRIQS/triqs/pull/799
#include <ostream>
#include <sstream>

#include <triqs/utility/variant_extensions.hpp>
#include <triqs/test_tools/arrays.hpp>

#include <realevol/time_expr.hpp>

#include <realevol/operators/many_body_operator.hpp>
#include <realevol/hilbert_space/hilbert_space.hpp>
#include <realevol/hilbert_space/imperative_operator.hpp>
#include <realevol/hilbert_space/state.hpp>
#include <realevol/hilbert_space/state_view.hpp>

using namespace realevol;
namespace hs = realevol::hilbert_space;

template<typename T> std::string as_string(T x) {
 std::stringstream ss; ss << x;
 return ss.str();
}

template<typename StateType>
void check_state(StateType const& st, std::map<int,typename StateType::value_type> const& ref) {
 foreach(st, [&st,&ref](int i, typename StateType::value_type a){
  auto it = ref.find(st.get_hilbert().get_fock_state(i));
  if(it != ref.end()) { EXPECT_CLOSE(it->second, a); }
  else                { EXPECT_EQ(.0, a); }
 });
}

TEST(hilbert_space, fundamental_operator_set) {
 using hs::statistic_enum::Fermion;
 using hs::statistic_enum::Boson;

 hs::fundamental_operator_set fop1(std::vector<int>(2,4));
 for (int i=0; i<2; ++i)
  for (int j=0; j<4; ++j) {
   EXPECT_EQ(4*i + j, (fop1[{i,j}]));
  }
 EXPECT_EQ(8, fop1.size());

 hs::fundamental_operator_set fop2;
 fop2 = fop1;
 EXPECT_EQ(5, (fop2[{1,1}]));

 hs::fundamental_operator_set fop3({},std::set<int>{1,2,3});
 for (int i=0; i<4; ++i) fop3.insert_fermion(i);
 EXPECT_EQ(2, (fop3[{2}]));
 EXPECT_EQ(3, fop3.pos({3}, Fermion));
 EXPECT_EQ(1, fop3.pos({2}, Boson));
 EXPECT_EQ(4, fop3.size(Fermion));
 EXPECT_EQ(3, fop3.size(Boson));
 EXPECT_EQ(7, fop3.size());

 hs::fundamental_operator_set fop4;
 for (int i=0; i<2; ++i) fop4.insert_fermion("up",i);
 for (int i=0; i<2; ++i) fop4.insert_fermion("down",i);
 EXPECT_EQ(0, (fop4[{"down",0}]));
 EXPECT_EQ(4, fop4.size());

 for (int i=0; i<2; ++i) fop4.insert_boson("B",i);
 EXPECT_EQ(0, fop4.pos({"down",0}, Fermion));
 EXPECT_EQ(1, fop4.pos({"B",1}, Boson));
 EXPECT_EQ(4, fop4.size(Fermion));
 EXPECT_EQ(2, fop4.size(Boson));
 EXPECT_EQ(6, fop4.size());

 h5::file fops_file("fops.h5", 'w');
 h5_write_attribute(fops_file, "fop", fop4);
 hs::fundamental_operator_set fop5;
 h5_read_attribute(fops_file, "fop", fop5);
 EXPECT_EQ(fop4, fop5);
}

TEST(hilbert_space, hilbert_space) {
 hs::fundamental_operator_set fop(std::vector<int>(2,4),std::vector<int>(1,2));

 hs::hilbert_space hs1(fop, {2,2});
 EXPECT_EQ(4096, hs1.size());

 EXPECT_TRUE(hs1.has_state(130));        // Fock state 130 is here
 EXPECT_FALSE(hs1.has_state(4096));      // Fock state 4096 is not here
 EXPECT_EQ(520, hs1.get_fock_state(520)); // Fock state for index 520
 EXPECT_EQ(520, hs1.get_state_index(hs1.get_fock_state(520))); // index of fock state 120

 std::set<hs::indices_t> ind_f = {};
 std::multiset<hs::indices_t> ind_b = {};
 // Fock state for vacuum
 EXPECT_EQ(0, hs1.get_fock_state(fop, ind_f));
 EXPECT_EQ(0, hs1.get_fock_state(fop, ind_f, ind_b));
 // Fock state for c^+(0,1)c^+(1,3)|vac>
 ind_f = {{0,1},{1,3}};
 EXPECT_EQ(130, hs1.get_fock_state(fop,ind_f));
 // Fock state for c^+(0,1)c^+(1,3)a^+(0,0)a^+(0,1)a^+(0,1)|vac>
 ind_b = {{0,0}, {0,1}, {0,1}};
 EXPECT_EQ(2434, hs1.get_fock_state(fop,ind_f,ind_b));

 hs::hilbert_space hs2;
 hs2 = hs1;
 EXPECT_EQ(4096, hs2.size());

 // HDF5
 auto hs_h5 = rw_h5(hs1, "hilbert_space");
 EXPECT_EQ(hs1, hs_h5);
}

TEST(hilbert_space, fock_state) {
 hs::fundamental_operator_set fop;
 for (int i=0; i<4; ++i) fop.insert_fermion(i);
 for (int i=0; i<3; ++i) fop.insert_boson(i);

 hs::hilbert_space hs(fop, {1,1,1});

 hs::fock_state_t fs1 = hs.get_fock_state(10);
 EXPECT_EQ(10, fs1);
 hs::fock_state_t fs2 = fs1;
 EXPECT_EQ(10, fs2);
}

TEST(hilbert_space, state) {
 hs::fundamental_operator_set fop;
 for (int i=0; i<5; ++i) fop.insert_fermion("up",i);
 for (int i=0; i<3; ++i) fop.insert_boson("B",i);

 hs::hilbert_space h_full(fop, {1,1,1});

 hs::state<hs::hilbert_space,double, true> st(h_full);
 st(0) = 3.0;
 st(3) = 5.0;
 st(77) = 1.0;
 check_state(st, {{3,5.0},{0,3.0},{77,1}});
}

TEST(hilbert_space, state_view) {
 hs::fundamental_operator_set fop;
 for (int i=0; i<5; ++i) fop.insert_fermion("up",i);
 for (int i=0; i<3; ++i) fop.insert_boson("B",i);

 hs::hilbert_space h_full(fop, {1,1,1});
 triqs::arrays::vector<double> amplitudes(h_full.size());
 amplitudes() = 0;

 hs::state_view<hs::hilbert_space,double> sv(amplitudes, h_full);
 sv(0) = 3.0;
 sv(3) = 5.0;
 sv(77) = 1.0;
 check_state(sv, {{0,3.0},{3,5.0},{77,1}});
 amplitudes[2] = 4.0;
 check_state(sv, {{0,3.0},{2,4.0},{3,5.0},{77,1}});
}

TEST(hilbert_space, imperative_operator) {
 hs::fundamental_operator_set fop;
 for (int i=0; i<4; ++i) fop.insert_fermion("up",i);
 fop.insert_boson("B",0);
 fop.insert_boson("B",1);

 hs::hilbert_space h_full(fop, {2,2});

 using operators::c;
 using operators::c_dag;
 using operators::a;
 using operators::a_dag;

 auto H = 3 * c_dag("up",1) * c("up",1) + 2 * c_dag("up",2) * c("up",2) + c("up",1) * c("up",2);
 EXPECT_EQ("3*C^+(up,1)C(up,1) + 2*C^+(up,2)C(up,2) + -1*C(up,2)C(up,1)", as_string(H));

 auto opH = hs::imperative_operator<hs::hilbert_space>(H, fop, h_full);

 hs::state<hs::hilbert_space, double, true> old_state(h_full);
 old_state(7) = 1.0;
 check_state(old_state, {{7,1.0}});

 auto new_state = opH(old_state);
 check_state(new_state, {{1,-1.0}, {7,5.0}});

 auto Hb = 3.0 * c_dag("up",0) * c("up",1) * a_dag("B",0) * a_dag("B",0) * a("B",1) +
           3.0 * c_dag("up",1) * c("up",3) * a_dag("B",0) * a("B",1) * a("B",1);

 EXPECT_EQ("3*C^+(up,0)C(up,1)[A^+(B,0)]^2A(B,1) + 3*C^+(up,1)C(up,3)A^+(B,0)[A(B,1)]^2", as_string(Hb));

 old_state(7) = 0;
 old_state(157) = 0.5; // |\up0>|\up2>|\up3>|1>_B0|2>_B1
 old_state(206) = 0.5; // |\up1>|\up2>|\up3>|0>_B0|3>_B1
 check_state(old_state, {{157,0.5}, {206,0.5}});

 auto opHb = hs::imperative_operator<hs::hilbert_space>(Hb, fop, h_full);

 new_state = opHb(old_state);
 check_state(new_state, {{39,-3.0}, {173,1.5*std::sqrt(6.0)}});
}

TEST(hilbert_space, sub_hilbert_space) {
 using operators::c_dag;
 using operators::a;
 auto CdagA = c_dag("up",0) * a("B",1);
 EXPECT_EQ("1*C^+(up,0)A(B,1)", as_string(CdagA));

 hs::fundamental_operator_set fop;
 for (int i=0; i<2; ++i) fop.insert_fermion("down",i);
 for (int i=0; i<2; ++i) fop.insert_fermion("up",i);
 for (int i=0; i<2; ++i) fop.insert_boson("B",i);

 hs::hilbert_space hs(fop,{2,2});

 hs::sub_hilbert_space phs0(0);
 phs0.add_fock_state(hs.get_fock_state(64)); // 01 00 0000
 phs0.add_fock_state(hs.get_fock_state(65)); // 01 00 0001
 phs0.add_fock_state(hs.get_fock_state(66)); // 01 00 0010
 phs0.add_fock_state(hs.get_fock_state(67)); // 01 00 0011

 EXPECT_TRUE(phs0.has_state(hs.get_fock_state(66))); // phs0 has state 66
 EXPECT_FALSE(phs0.has_state(hs.get_fock_state(6))); // phs0 has state 6

 hs::sub_hilbert_space phs1(1);
 phs1.add_fock_state(hs.get_fock_state(4)); // 00 00 0100
 phs1.add_fock_state(hs.get_fock_state(5)); // 00 00 0101
 phs1.add_fock_state(hs.get_fock_state(6)); // 00 00 0110
 phs1.add_fock_state(hs.get_fock_state(7)); // 00 00 0111

 EXPECT_FALSE(phs1.has_state(hs.get_fock_state(66))); // phs1 has state 66
 EXPECT_TRUE(phs1.has_state(hs.get_fock_state(6)));   // phs1 has state 6

 std::vector<int> Cdagmap(2,-1);
 Cdagmap[phs0.get_index()] = phs1.get_index();
 std::vector<hs::sub_hilbert_space> sub1{phs0, phs1};
 auto opCdagA = hs::imperative_operator<hs::sub_hilbert_space, double, true>(CdagA, fop, hs, Cdagmap, &sub1);

 hs::state<hs::sub_hilbert_space,double,false> start(phs0);
 start(0) = 1.0;
 start(1) = 2.0;
 start(2) = 3.0;
 start(3) = 4.0;

 check_state(start, {{64,1.0}, {65,2.0}, {66,3.0}, {67,4.0}});
 check_state(opCdagA(start), {{4,1.0}, {5,-2.0}, {6,-3.0}, {7,4.0}});

 // HDF5
 auto hs_h5 = rw_h5(phs1, "sub_hilbert_space");
 EXPECT_EQ(phs1, hs_h5);
}

TEST(hilbert_space, QuarticOperators) {
 hs::fundamental_operator_set fops;
 fops.insert_fermion("down",0);
 fops.insert_fermion("up",0);
 fops.insert_fermion("down",1);
 fops.insert_fermion("up",1);
 fops.insert_boson("B",0);
 fops.insert_boson("B",1);

 hs::hilbert_space hs(fops,{2,2});
 EXPECT_EQ(256, hs.size());

 using operators::c;
 using operators::c_dag;
 using operators::a;
 using operators::a_dag;

 auto quartic_op1 = -1.0*c_dag("up",0)*c_dag("down",1)*c("up",1)*c("down",0);

 hs::state<hs::hilbert_space,double, false> st1(hs);
 st1(9) = 1.0; // 00 00 1001
 check_state(st1, {{9,1.0}}); // old state
 EXPECT_EQ("1*C^+(down,1)C^+(up,0)C(up,1)C(down,0)", as_string(quartic_op1)); // quartic operator 1
 check_state(hs::imperative_operator<hs::hilbert_space>(quartic_op1,fops,hs)(st1),
             {{6,1.0}}/* 00 00 0110 */); // new state

 auto quartic_op2 = 1.0*c_dag("up",0)*c("down",1)*a_dag("B",0)*a("B",1);
 hs::state<hs::hilbert_space,double, false> st2(hs);
 st2(218) = 1.0; // 11 01 1010
 check_state(st2, {{218,1.0}}); // old state
 EXPECT_EQ("1*C^+(up,0)C(down,1)A^+(B,0)A(B,1)", as_string(quartic_op2)); // quartic operator 2
 check_state(hs::imperative_operator<hs::hilbert_space>(quartic_op2,fops,hs)(st2),
             {{172,std::sqrt(6.0)}}/* 10 10 1100 */); // new state
}

TEST(hilbert_space, StateProjection) {
 hs::fundamental_operator_set fop;
 for (int i=0; i<3; ++i) fop.insert_fermion("s",i);
 for (int i=0; i<3; ++i) fop.insert_boson("b",i);

 hs::hilbert_space hs_full(fop,{2,2,2});

 hs::state<hs::hilbert_space,double,true> st(hs_full);
 st(0) = 0.1;
 st(2) = 0.2;
 st(4) = 0.3;
 st(6) = 0.4;
 st(12) = 0.5;
 check_state(st, {{12,0.5},{0,0.1},{2,0.2},{4,0.3},{6,0.4}}); // original state

 hs::sub_hilbert_space hs(0);
 hs.add_fock_state(hs_full.get_fock_state(4));
 hs.add_fock_state(hs_full.get_fock_state(5));
 hs.add_fock_state(hs_full.get_fock_state(6));
 hs.add_fock_state(hs_full.get_fock_state(7));
 hs.add_fock_state(hs_full.get_fock_state(14));

 auto proj_st = hs::project<hs::state<hs::sub_hilbert_space,double,false>>(st,hs);
 check_state (proj_st, {{4,0.3},{6,0.4}}); // projected state
}

TEST(hilbert_space, time_expr_real) {
 hs::fundamental_operator_set fops;
 fops.insert_fermion("up",0);
 fops.insert_fermion("down",0);
 fops.insert_fermion("up",1);
 fops.insert_fermion("down",1);
 fops.insert_boson("B",0);
 fops.insert_boson("B",1);

 hs::hilbert_space hs(fops,{2,2});
 EXPECT_EQ(256, hs.size());

 operators::many_body_operator_generic<time_expr> op;

 using operators::c;
 using operators::c_dag;
 using operators::a;
 using operators::a_dag;

 op = 1.0*c_dag<time_expr>("up",0)*c<time_expr>("down",1)
         *a_dag<time_expr>("B",0)*a<time_expr>("B",1);

 // Spin-flips
 op += "2*t^2"*c_dag<time_expr>("down",0)*c<time_expr>("up",1);
 op += "2*t^2"*c_dag<time_expr>("up",0)*c<time_expr>("down",1);

 double t = 0.2;

 hs::state<hs::hilbert_space,dcomplex,false> st1(hs);
 st1(218) = 1.0; // 11 01 1010
 check_state(st1, {{218,1.0}});
 EXPECT_EQ("2*t^2*C^+(down,0)C(up,1) + 2*t^2*C^+(up,0)C(down,1) + 1*C^+(up,0)C(down,1)A^+(B,0)A(B,1)",
           as_string(op));
 auto imp_op = hs::imperative_operator<hs::hilbert_space,time_expr>(op,fops,hs);

 check_state(imp_op(st1,t), {{172,std::sqrt(6.0)},{211,-0.08},{220,0.08}});

 // precompute imp_op at t=0.3
 t = 0.3;
 imp_op.update_coeffs([t](time_expr & M){ M = M(t);});
 check_state(imp_op(st1,t), {{172,std::sqrt(6.0)},{211,-0.18},{220,0.18}});
}

TEST(hilbert_space, time_expr_complex) {
 hs::fundamental_operator_set fops;
 fops.insert_fermion("up",0);
 fops.insert_fermion("down",0);
 fops.insert_fermion("up",1);
 fops.insert_fermion("down",1);
 fops.insert_boson("B",0);
 fops.insert_boson("B",1);

 hs::hilbert_space hs(fops,{2,2});
 EXPECT_EQ(256, hs.size());

 operators::many_body_operator_generic<time_expr> op;

 using operators::c;
 using operators::c_dag;
 using operators::a;
 using operators::a_dag;

 op = 1.0*c_dag<time_expr>("up",0)*c<time_expr>("down",1)
         *a_dag<time_expr>("B",0)*a<time_expr>("B",1);

 // Spin-flips
 op += 1i*("2*t^2"*c_dag<time_expr>("down",0)*c<time_expr>("up",1));
 op += 1i*("2*t^2"*c_dag<time_expr>("up",0)*c<time_expr>("down",1));

 double t = 0.2;

 hs::state<hs::hilbert_space,dcomplex,false> st1(hs);
 st1(218) = 1.0; // 11 01 1010
 check_state(st1, {{218, 1.0}});
 EXPECT_EQ("(0,2*t^2)*C^+(down,0)C(up,1) + (0,2*t^2)*C^+(up,0)C(down,1) + 1*C^+(up,0)C(down,1)A^+(B,0)A(B,1)",
           as_string(op));
 auto imp_op = hs::imperative_operator<hs::hilbert_space,time_expr>(op,fops,hs);

 check_state(imp_op(st1,t), {{172, std::sqrt(6.0)},{211,-0.08i},{220,0.08i}});

 // precompute imp_op at t=0.3
 t = 0.3;
 imp_op.update_coeffs([t](time_expr & M){ M = M(t);});
 check_state(imp_op(st1,t), {{172, std::sqrt(6.0)},{211,-0.18i},{220,0.18i}});
}

MAKE_MAIN;
