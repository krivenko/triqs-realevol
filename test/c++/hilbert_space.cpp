#include <triqs/test_tools/gfs.hpp>
#include <sstream>

#include <time_expr.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include <triqs/hilbert_space/state.hpp>

using namespace realevol;
using namespace triqs::hilbert_space;

#define EXPECT_PRINT(X, Y) {std::stringstream ss; ss << Y; EXPECT_EQ(X,ss.str()); }
#define ASSERT_PRINT(X, Y) {std::stringstream ss; ss << Y; ASSERT_EQ(X,ss.str()); }

auto I = dcomplex(0,1.0);

TEST(hilbert_space, fundamental_operator_set) {
 using triqs::hilbert_space::statistic_enum::Fermion;
 using triqs::hilbert_space::statistic_enum::Boson;

 fundamental_operator_set fop1(std::vector<int>(2,4));
 for (int i=0; i<2; ++i)
  for (int j=0; j<4; ++j) {
   EXPECT_EQ(4*i + j, (fop1[{i,j}]));
  }
 EXPECT_EQ(8, fop1.size());

 fundamental_operator_set fop2;
 fop2 = fop1;
 EXPECT_EQ(5, (fop2[{1,1}]));

 fundamental_operator_set fop3({},std::set<int>{1,2,3});
 for (int i=0; i<4; ++i) fop3.insert_fermion(i);
 EXPECT_EQ(2, (fop3[{2}]));
 EXPECT_EQ(3, fop3.pos({3}, Fermion));
 EXPECT_EQ(1, fop3.pos({2}, Boson));
 EXPECT_EQ(4, fop3.size(Fermion));
 EXPECT_EQ(3, fop3.size(Boson));
 EXPECT_EQ(7, fop3.size());

 fundamental_operator_set fop4;
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
}

TEST(hilbert_space, hilbert_space) {
 fundamental_operator_set fop(std::vector<int>(2,4),std::vector<int>(1,2));

 using triqs::hilbert_space::hilbert_space;
 hilbert_space hs1(fop, {2,2});
 EXPECT_EQ(4096, hs1.size());

 EXPECT_TRUE(hs1.has_state(130));        // Fock state 130 is here
 EXPECT_FALSE(hs1.has_state(4096));      // Fock state 4096 is not here
 EXPECT_EQ(520, hs1.get_fock_state(520)); // Fock state for index 520
 EXPECT_EQ(520, hs1.get_state_index(hs1.get_fock_state(520))); // index of fock state 120

 std::set<fundamental_operator_set::indices_t> ind_f = {};
 std::multiset<fundamental_operator_set::indices_t> ind_b = {};
 // Fock state for vacuum
 EXPECT_EQ(0, hs1.get_fock_state(fop, ind_f));
 EXPECT_EQ(0, hs1.get_fock_state(fop, ind_f, ind_b));
 // Fock state for c^+(0,1)c^+(1,3)|vac>
 ind_f = {{0,1},{1,3}};
 EXPECT_EQ(130, hs1.get_fock_state(fop,ind_f));
 // Fock state for c^+(0,1)c^+(1,3)a^+(0,0)a^+(0,1)a^+(0,1)|vac>
 ind_b = {{0,0}, {0,1}, {0,1}};
 EXPECT_EQ(2434, hs1.get_fock_state(fop,ind_f,ind_b));

 hilbert_space hs2;
 hs2 = hs1;
 EXPECT_EQ(4096, hs2.size());

 // HDF5
 auto hs_h5 = rw_h5(hs1, "hilbert_space");
 EXPECT_EQ(hs1, hs_h5);
}

TEST(hilbert_space, fock_state) {
 fundamental_operator_set fop;
 for (int i=0; i<4; ++i) fop.insert_fermion(i);
 for (int i=0; i<3; ++i) fop.insert_boson(i);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space hs(fop, {1,1,1});

 fock_state_t fs1 = hs.get_fock_state(10);
 EXPECT_EQ(10, fs1);
 fock_state_t fs2 = fs1;
 EXPECT_EQ(10, fs2);
}

TEST(hilbert_space, state) {
 fundamental_operator_set fop;
 for (int i=0; i<5; ++i) fop.insert_fermion("up",i);
 for (int i=0; i<3; ++i) fop.insert_boson("B",i);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space h_full(fop, {1,1,1});

 state<hilbert_space,double, true> st(h_full);
 st(0) = 3.0;
 st(3) = 5.0;
 st(77) = 1.0;
 EXPECT_PRINT(" +(1)|77> +(5)|3> +(3)|0>", st);
}

TEST(hilbert_space, imperative_operator) {
 fundamental_operator_set fop;
 for (int i=0; i<4; ++i) fop.insert_fermion("up",i);
 fop.insert_boson("B",0);
 fop.insert_boson("B",1);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space h_full(fop, {2,2});

 using triqs::operators::c;
 using triqs::operators::c_dag;
 using triqs::operators::a;
 using triqs::operators::a_dag;

 auto H = 3 * c_dag("up",1) * c("up",1) + 2 * c_dag("up",2) * c("up",2) + c("up",1) * c("up",2);
 EXPECT_PRINT("3*C^+(up,1)C(up,1) + 2*C^+(up,2)C(up,2) + -1*C(up,2)C(up,1)", H);

 auto opH = imperative_operator<hilbert_space>(H, fop, h_full);

 state<hilbert_space, double, true> old_state(h_full);
 old_state(7) = 1.0;
 EXPECT_PRINT(" +(1)|7>", old_state);

 auto new_state = opH(old_state);
 EXPECT_PRINT(" +(-1)|1> +(5)|7>",new_state);

 auto Hb = 3.0 * c_dag("up",0) * c("up",1) * a_dag("B",0) * a_dag("B",0) * a("B",1) +
           3.0 * c_dag("up",1) * c("up",3) * a_dag("B",0) * a("B",1) * a("B",1);

 EXPECT_PRINT("3*C^+(up,0)C(up,1)[A^+(B,0)]^2A(B,1) + 3*C^+(up,1)C(up,3)A^+(B,0)[A(B,1)]^2", Hb);

 old_state(7) = 0;
 old_state(157) = 0.5; // |\up0>|\up2>|\up3>|1>_B0|2>_B1
 old_state(206) = 0.5; // |\up1>|\up2>|\up3>|0>_B0|3>_B1
 EXPECT_PRINT(" +(0.5)|206> +(0.5)|157>", old_state);

 auto opHb = imperative_operator<hilbert_space>(Hb, fop, h_full);

 new_state = opHb(old_state);
 EXPECT_PRINT(" +(-3)|39> +(3.67423)|173>", new_state);
}

TEST(hilbert_space, sub_hilbert_space) {
 using triqs::operators::c_dag;
 using triqs::operators::a;
 auto CdagA = c_dag("up",0) * a("B",1);
 EXPECT_PRINT("1*C^+(up,0)A(B,1)", CdagA);

 fundamental_operator_set fop;
 for (int i=0; i<2; ++i) fop.insert_fermion("down",i);
 for (int i=0; i<2; ++i) fop.insert_fermion("up",i);
 for (int i=0; i<2; ++i) fop.insert_boson("B",i);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space hs(fop,{2,2});

 sub_hilbert_space phs0(0);
 phs0.add_fock_state(hs.get_fock_state(64)); // 01 00 0000
 phs0.add_fock_state(hs.get_fock_state(65)); // 01 00 0001
 phs0.add_fock_state(hs.get_fock_state(66)); // 01 00 0010
 phs0.add_fock_state(hs.get_fock_state(67)); // 01 00 0011

 EXPECT_TRUE(phs0.has_state(hs.get_fock_state(66))); // phs0 has state 66
 EXPECT_FALSE(phs0.has_state(hs.get_fock_state(6))); // phs0 has state 6

 sub_hilbert_space phs1(1);
 phs1.add_fock_state(hs.get_fock_state(4)); // 00 00 0100
 phs1.add_fock_state(hs.get_fock_state(5)); // 00 00 0101
 phs1.add_fock_state(hs.get_fock_state(6)); // 00 00 0110
 phs1.add_fock_state(hs.get_fock_state(7)); // 00 00 0111

 EXPECT_FALSE(phs1.has_state(hs.get_fock_state(66))); // phs1 has state 66
 EXPECT_TRUE(phs1.has_state(hs.get_fock_state(6)));   // phs1 has state 6

 std::vector<int> Cdagmap(2,-1);
 Cdagmap[phs0.get_index()] = phs1.get_index();
 std::vector<sub_hilbert_space> sub1{phs0, phs1};
 auto opCdagA = imperative_operator<sub_hilbert_space, double, true>(CdagA, fop, hs, Cdagmap, &sub1);

 state<sub_hilbert_space,double,false> start(phs0);
 start(0) = 1.0;
 start(1) = 2.0;
 start(2) = 3.0;
 start(3) = 4.0;

 EXPECT_PRINT(" +(1)|64> +(2)|65> +(3)|66> +(4)|67>", start);
 EXPECT_PRINT(" +(1)|4> +(-2)|5> +(-3)|6> +(4)|7>", opCdagA(start));

 // HDF5
 auto hs_h5 = rw_h5(phs1, "sub_hilbert_space");
 EXPECT_EQ(phs1, hs_h5);
}

TEST(hilbert_space, QuarticOperators) {
 fundamental_operator_set fops;
 fops.insert_fermion("down",0);
 fops.insert_fermion("up",0);
 fops.insert_fermion("down",1);
 fops.insert_fermion("up",1);
 fops.insert_boson("B",0);
 fops.insert_boson("B",1);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space hs(fops,{2,2});
 EXPECT_EQ(256, hs.size());

 using triqs::operators::c;
 using triqs::operators::c_dag;
 using triqs::operators::a;
 using triqs::operators::a_dag;

 auto quartic_op1 = -1.0*c_dag("up",0)*c_dag("down",1)*c("up",1)*c("down",0);

 state<hilbert_space,double, false> st1(hs);
 st1(9) = 1.0; // 00 00 1001
 EXPECT_PRINT(" +(1)|9>", st1); // old state
 EXPECT_PRINT("1*C^+(down,1)C^+(up,0)C(up,1)C(down,0)", quartic_op1); // quartic operator 1
 EXPECT_PRINT(" +(1)|6>"/* 00 00 0110 */,
              imperative_operator<hilbert_space>(quartic_op1,fops,hs)(st1)); // new state

 auto quartic_op2 = 1.0*c_dag("up",0)*c("down",1)*a_dag("B",0)*a("B",1);
 state<hilbert_space,double, false> st2(hs);
 st2(218) = 1.0; // 11 01 1010
 EXPECT_PRINT(" +(1)|218>", st2); // old state
 EXPECT_PRINT("1*C^+(up,0)C(down,1)A^+(B,0)A(B,1)", quartic_op2); // quartic operator 2
 EXPECT_PRINT(" +(2.44949)|172>"/* 10 10 1100 */,
              imperative_operator<hilbert_space>(quartic_op2,fops,hs)(st2)); // new state
}

TEST(hilbert_space, StateProjection) {
 fundamental_operator_set fop;
 for (int i=0; i<3; ++i) fop.insert_fermion("s",i);
 for (int i=0; i<3; ++i) fop.insert_boson("b",i);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space hs_full(fop,{2,2,2});

 state<hilbert_space,double,true> st(hs_full);
 st(0) = 0.1;
 st(2) = 0.2;
 st(4) = 0.3;
 st(6) = 0.4;
 st(12) = 0.5;
 EXPECT_PRINT(" +(0.5)|12> +(0.1)|0> +(0.2)|2> +(0.3)|4> +(0.4)|6>", st); // original state

 sub_hilbert_space hs(0);
 hs.add_fock_state(hs_full.get_fock_state(4));
 hs.add_fock_state(hs_full.get_fock_state(5));
 hs.add_fock_state(hs_full.get_fock_state(6));
 hs.add_fock_state(hs_full.get_fock_state(7));
 hs.add_fock_state(hs_full.get_fock_state(14));

 auto proj_st = project<state<sub_hilbert_space,double,false>>(st,hs);
 EXPECT_PRINT(" +(0.3)|4> +(0.4)|6>", proj_st); // projected state
}

TEST(hilbert_space, time_expr_real) {
 fundamental_operator_set fops;
 fops.insert_fermion("up",0);
 fops.insert_fermion("down",0);
 fops.insert_fermion("up",1);
 fops.insert_fermion("down",1);
 fops.insert_boson("B",0);
 fops.insert_boson("B",1);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space hs(fops,{2,2});
 EXPECT_EQ(256, hs.size());

 triqs::operators::many_body_operator_generic<time_expr> op;

 using triqs::operators::c;
 using triqs::operators::c_dag;
 using triqs::operators::a;
 using triqs::operators::a_dag;

 op = 1.0*c_dag<time_expr>("up",0)*c<time_expr>("down",1)
         *a_dag<time_expr>("B",0)*a<time_expr>("B",1);

 // Spin-flips
 op += "2*t^2"*c_dag<time_expr>("down",0)*c<time_expr>("up",1);
 op += "2*t^2"*c_dag<time_expr>("up",0)*c<time_expr>("down",1);

 double t = 0.2;

 state<hilbert_space,dcomplex,false> st1(hs);
 st1(218) = 1.0; // 11 01 1010
 EXPECT_PRINT(" +((1,0))|218>", st1);
 EXPECT_PRINT("2*t^2*C^+(down,0)C(up,1) + 2*t^2*C^+(up,0)C(down,1) + 1*C^+(up,0)C(down,1)A^+(B,0)A(B,1)", op);
 auto imp_op = imperative_operator<hilbert_space,time_expr>(op,fops,hs);

 EXPECT_PRINT(" +((2.44949,0))|172> +((-0.08,0))|211> +((0.08,0))|220>", imp_op(st1,t));

 // precompute imp_op at t=0.3
 t = 0.3;
 imp_op.update_coeffs([t](time_expr & M){ M = M(t);});
 EXPECT_PRINT(" +((2.44949,0))|172> +((-0.18,0))|211> +((0.18,0))|220>", imp_op(st1,t));
}

TEST(hilbert_space, time_expr_complex) {
 fundamental_operator_set fops;
 fops.insert_fermion("up",0);
 fops.insert_fermion("down",0);
 fops.insert_fermion("up",1);
 fops.insert_fermion("down",1);
 fops.insert_boson("B",0);
 fops.insert_boson("B",1);

 using triqs::hilbert_space::hilbert_space;
 hilbert_space hs(fops,{2,2});
 EXPECT_EQ(256, hs.size());

 triqs::operators::many_body_operator_generic<time_expr> op;

 using triqs::operators::c;
 using triqs::operators::c_dag;
 using triqs::operators::a;
 using triqs::operators::a_dag;

 op = 1.0*c_dag<time_expr>("up",0)*c<time_expr>("down",1)
         *a_dag<time_expr>("B",0)*a<time_expr>("B",1);

 // Spin-flips
 op += I*("2*t^2"*c_dag<time_expr>("down",0)*c<time_expr>("up",1));
 op += I*("2*t^2"*c_dag<time_expr>("up",0)*c<time_expr>("down",1));

 double t = 0.2;

 state<hilbert_space,dcomplex,false> st1(hs);
 st1(218) = 1.0; // 11 01 1010
 EXPECT_PRINT(" +((1,0))|218>", st1);
 EXPECT_PRINT("(0,2*t^2)*C^+(down,0)C(up,1) + (0,2*t^2)*C^+(up,0)C(down,1) + 1*C^+(up,0)C(down,1)A^+(B,0)A(B,1)", op);
 auto imp_op = imperative_operator<hilbert_space,time_expr>(op,fops,hs);

 EXPECT_PRINT(" +((2.44949,0))|172> +((0,-0.08))|211> +((0,0.08))|220>", imp_op(st1,t));

 // precompute imp_op at t=0.3
 t = 0.3;
 imp_op.update_coeffs([t](time_expr & M){ M = M(t);});
 EXPECT_PRINT(" +((2.44949,0))|172> +((0,-0.18))|211> +((0,0.18))|220>", imp_op(st1,t));
}

MAKE_MAIN;
