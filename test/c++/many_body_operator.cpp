#include <triqs/utility/first_include.hpp>
#include <triqs/test_tools/arrays.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/complex.hpp>

#include <time_expr.hpp>
#include <triqs/operators/many_body_operator.hpp>

using triqs::operators::many_body_operator;
using triqs::operators::c_dag;
using triqs::operators::c;
using triqs::operators::n;
using triqs::operators::a_dag;
using triqs::operators::a;

using namespace realevol;
using std::to_string;

#define EXPECT_PRINT(X, Y) {std::stringstream ss; ss << Y; EXPECT_EQ(X,ss.str()); }
#define ASSERT_PRINT(X, Y) {std::stringstream ss; ss << Y; ASSERT_EQ(X,ss.str()); }

auto I = dcomplex(0,1.0);

TEST(operator, real_or_complex) {

 // Operators without indices
 auto op_with_no_indices = c() + c_dag() - n();
 EXPECT_PRINT("1*C^+() + 1*C() + -1*C^+()C()", op_with_no_indices);

 // Operators with many indices
 auto op_with_many_indices = c(1,20,"a",true,-2) + c_dag(3,15,"b",false,-5);
 EXPECT_PRINT("1*C^+(3,15,b,0,-5) + 1*C(1,20,a,1,-2)", op_with_many_indices);

 // Commutation relations
 std::string ref;
 for(int i = 0; i < 3; ++i)
 for(int j = 0; j < 3; ++j) {
  ref = (i == j ? "1" : "0");
  EXPECT_PRINT(ref, c_dag(i)*c(j) + c(j)*c_dag(i));
  ref = std::string(i == j ? "-1 + " : "") + "2*C^+(" + to_string(i) + ")C(" + to_string(j) + ")";
  EXPECT_PRINT(ref, c_dag(i)*c(j) - c(j)*c_dag(i));
  ref = (i == j ? "1" : "0");
  EXPECT_PRINT(ref, a(i)*a_dag(j) - a_dag(j)*a(i));
  ref = std::string(i == j ? "1 + " : "") + "2*A^+(" + to_string(j) + ")A(" + to_string(i) + ")";
  EXPECT_PRINT(ref, a(i)*a_dag(j) + a_dag(j)*a(i));
 }

 // Algebra
 auto C = c(0), Cd = c_dag(1);
 auto A = a(0), Ad = a_dag(1);

 EXPECT_PRINT("1*C(0)",   C);
 EXPECT_PRINT("1*C^+(1)", Cd);
 EXPECT_PRINT("1*A(0)",   A);
 EXPECT_PRINT("1*A^+(1)", Ad);

 // Unary minus
 EXPECT_PRINT("-1*C(0)",   -C);
 EXPECT_PRINT("-1*C^+(1)", -Cd);
 EXPECT_PRINT("-1*A(0)",   -A);
 EXPECT_PRINT("-1*A^+(1)", -Ad);

 // Addition
 EXPECT_PRINT("2 + 1*C(0)",       C + 2.0);
 EXPECT_PRINT("2 + 1*C^+(1)",     Cd + 2.0);
 EXPECT_PRINT("2 + 1*A(0)",       A + 2.0);
 EXPECT_PRINT("2 + 1*A^+(1)",     Ad + 2.0);
 EXPECT_PRINT("(0,2) + 1*C(0)",   C + 2.0*I);
 EXPECT_PRINT("(0,2) + 1*C^+(1)", Cd + 2.0*I);
 EXPECT_PRINT("(0,2) + 1*A(0)",   A + 2.0*I);
 EXPECT_PRINT("(0,2) + 1*A^+(1)", Ad + 2.0*I);

 EXPECT_PRINT("2 + 1*C(0)",       2.0 + C);
 EXPECT_PRINT("2 + 1*C^+(1)",     2.0 + Cd);
 EXPECT_PRINT("2 + 1*A(0)",       2.0 + A);
 EXPECT_PRINT("2 + 1*A^+(1)",     2.0 + Ad);
 EXPECT_PRINT("(0,2) + 1*C(0)",   2.0*I + C);
 EXPECT_PRINT("(0,2) + 1*C^+(1)", 2.0*I + Cd);
 EXPECT_PRINT("(0,2) + 1*A(0)",   2.0*I + A);
 EXPECT_PRINT("(0,2) + 1*A^+(1)", 2.0*I + Ad);

 EXPECT_PRINT("1*C^+(1) + 1*C(0) + 1*A^+(1) + 1*A(0)", C + Cd + A + Ad);

 // Subtraction
 EXPECT_PRINT("-2 + 1*C(0)",       C - 2.0);
 EXPECT_PRINT("-2 + 1*C^+(1)",     Cd - 2.0);
 EXPECT_PRINT("-2 + 1*A(0)",       A - 2.0);
 EXPECT_PRINT("-2 + 1*A^+(1)",     Ad - 2.0);
 EXPECT_PRINT("(-0,-2) + 1*C(0)",   C - 2.0*I);
 EXPECT_PRINT("(-0,-2) + 1*C^+(1)", Cd - 2.0*I);
 EXPECT_PRINT("(-0,-2) + 1*A(0)",   A - 2.0*I);
 EXPECT_PRINT("(-0,-2) + 1*A^+(1)", Ad - 2.0*I);

 EXPECT_PRINT("2 + -1*C(0)",       2.0 - C);
 EXPECT_PRINT("2 + -1*C^+(1)",     2.0 - Cd);
 EXPECT_PRINT("2 + -1*A(0)",       2.0 - A);
 EXPECT_PRINT("2 + -1*A^+(1)",     2.0 - Ad);
 EXPECT_PRINT("(0,2) + -1*C(0)",   2.0*I - C);
 EXPECT_PRINT("(0,2) + -1*C^+(1)", 2.0*I - Cd);
 EXPECT_PRINT("(0,2) + -1*A(0)",   2.0*I - A);
 EXPECT_PRINT("(0,2) + -1*A^+(1)", 2.0*I - Ad);

 EXPECT_PRINT("-1*C^+(1) + 1*C(0) + -1*A^+(1) + -1*A(0)", C - Cd - A - Ad);

 // Multiplication
 EXPECT_PRINT("3*C(0)",       C * 3.0);
 EXPECT_PRINT("3*C^+(1)",     Cd * 3.0);
 EXPECT_PRINT("3*A(0)",       A * 3.0);
 EXPECT_PRINT("3*A^+(1)",     Ad * 3.0);
 EXPECT_PRINT("(0,3)*C(0)",   C * 3.0*I);
 EXPECT_PRINT("(0,3)*C^+(1)", Cd * 3.0*I);
 EXPECT_PRINT("(0,3)*A(0)",   A * 3.0*I);
 EXPECT_PRINT("(0,3)*A^+(1)", Ad * 3.0*I);

 EXPECT_PRINT("3*C(0)",       3.0 * C);
 EXPECT_PRINT("3*C^+(1)",     3.0 * Cd);
 EXPECT_PRINT("3*A(0)",       3.0 * A);
 EXPECT_PRINT("3*A^+(1)",     3.0 * Ad);
 EXPECT_PRINT("(0,3)*C(0)",   3.0*I * C);
 EXPECT_PRINT("(0,3)*C^+(1)", 3.0*I * Cd);
 EXPECT_PRINT("(0,3)*A(0)",   3.0*I * A);
 EXPECT_PRINT("(0,3)*A^+(1)", 3.0*I * Ad);

 EXPECT_PRINT("-2 + 2*C^+(2)C(2) + -2*C^+(2)A^+(2) + 2*C(2)A(2) + -1*[A^+(2)]^2 + 1*[A(2)]^2",
              (c(2) + c_dag(2) + a(2) + a_dag(2))*(c(2) - c_dag(2) + a(2) - a_dag(2)));

 // (n_up * a + n_dn * a^+)^3
 auto expr = n("up") * a(0) + n("dn") * a_dag(0);
 EXPECT_PRINT("1*C^+(dn)C(dn)A^+(0) + 1*C^+(up)C(up)A(0)", expr);
 expr = expr*expr*expr;
 EXPECT_PRINT("3*C^+(dn)C^+(up)C(up)C(dn)A^+(0) + 3*C^+(dn)C^+(up)C(up)C(dn)A(0) + "
              "1*C^+(dn)C(dn)[A^+(0)]^3 + 1*C^+(up)C(up)[A(0)]^3 + "
              "3*C^+(dn)C^+(up)C(up)C(dn)[A^+(0)]^2A(0) + 3*C^+(dn)C^+(up)C(up)C(dn)A^+(0)[A(0)]^2",
              expr);

 // Dagger
 auto X = I*c_dag(1) * c_dag(2) * c(3) * c(4) * a_dag(5) * a(6);
 EXPECT_PRINT("(0,-1)*C^+(1)C^+(2)C(4)C(3)A^+(5)A(6)", X);
 EXPECT_PRINT("(0,1)*C^+(3)C^+(4)C(2)C(1)A^+(6)A(5)", dagger(X));
}

TEST(operator, time_expr) {

 using te = realevol::time_expr;

 // Operators without indices
 auto op_with_no_indices = c<te>() + c_dag<te>() - n<te>();
 EXPECT_PRINT("1*C^+() + 1*C() + -1*C^+()C()", op_with_no_indices);

 // Operators with many indices
 auto op_with_many_indices = c<te>(1,20,"a",true,-2) + c_dag<te>(3,15,"b",false,-5);
 EXPECT_PRINT("1*C^+(3,15,b,0,-5) + 1*C(1,20,a,1,-2)", op_with_many_indices);

 // Commutation relations
 std::string ref;
 for(int i = 0; i < 3; ++i)
 for(int j = 0; j < 3; ++j) {
  ref = (i == j ? "1" : "0");
  EXPECT_PRINT(ref, c_dag<te>(i)*c<te>(j) + c<te>(j)*c_dag<te>(i));
  ref = std::string(i == j ? "-1 + " : "") + "2*C^+(" + to_string(i) + ")C(" + to_string(j) + ")";
  EXPECT_PRINT(ref, c_dag<te>(i)*c<te>(j) - c<te>(j)*c_dag<te>(i));
  ref = (i == j ? "1" : "0");
  EXPECT_PRINT(ref, a<te>(i)*a_dag<te>(j) - a_dag<te>(j)*a<te>(i));
  ref = std::string(i == j ? "1 + " : "") + "2*A^+(" + to_string(j) + ")A(" + to_string(i) + ")";
  EXPECT_PRINT(ref, a<te>(i)*a_dag<te>(j) + a_dag<te>(j)*a<te>(i));
 }

 // Algebra
 auto C = "t^2"_te * c<te>(0), Cd = c_dag<te>(1);
 auto A = a<te>(0), Ad = "sin(t)"_te * a_dag<te>(1);

 EXPECT_PRINT("t^2*C(0)",      C);
 EXPECT_PRINT("1*C^+(1)",      Cd);
 EXPECT_PRINT("1*A(0)",        A);
 EXPECT_PRINT("sin(t)*A^+(1)", Ad);

 // Unary minus
 EXPECT_PRINT("-(t^2)*C(0)",      -C);
 EXPECT_PRINT("-1*C^+(1)",        -Cd);
 EXPECT_PRINT("-1*A(0)",          -A);
 EXPECT_PRINT("-(sin(t))*A^+(1)", -Ad);

 // Addition
 EXPECT_PRINT("2 + t^2*C(0)",          C + 2.0);
 EXPECT_PRINT("2 + 1*C^+(1)",          Cd + 2.0);
 EXPECT_PRINT("2 + 1*A(0)",            A + 2.0);
 EXPECT_PRINT("2 + sin(t)*A^+(1)",     Ad + 2.0);
 EXPECT_PRINT("(0,2) + t^2*C(0)",      C + 2.0*I);
 EXPECT_PRINT("(0,2) + 1*C^+(1)",      Cd + 2.0*I);
 EXPECT_PRINT("(0,2) + 1*A(0)",        A + 2.0*I);
 EXPECT_PRINT("(0,2) + sin(t)*A^+(1)", Ad + 2.0*I);

 EXPECT_PRINT("2*t + t^2*C(0)",          C + "2*t"_te);
 EXPECT_PRINT("2*t + 1*C^+(1)",          Cd + "2*t"_te);
 EXPECT_PRINT("2*t + 1*A(0)",            A + "2*t"_te);
 EXPECT_PRINT("2*t + sin(t)*A^+(1)",     Ad + "2*t"_te);
 EXPECT_PRINT("(0,2*t) + t^2*C(0)",      C + I*"2*t"_te);
 EXPECT_PRINT("(0,2*t) + 1*C^+(1)",      Cd + I*"2*t"_te);
 EXPECT_PRINT("(0,2*t) + 1*A(0)",        A + I*"2*t"_te);
 EXPECT_PRINT("(0,2*t) + sin(t)*A^+(1)", Ad + I*"2*t"_te);

 EXPECT_PRINT("2 + t^2*C(0)",          2.0 + C);
 EXPECT_PRINT("2 + 1*C^+(1)",          2.0 + Cd);
 EXPECT_PRINT("2 + 1*A(0)",            2.0 + A);
 EXPECT_PRINT("2 + sin(t)*A^+(1)",     2.0 + Ad);
 EXPECT_PRINT("(0,2) + t^2*C(0)",      2.0*I + C);
 EXPECT_PRINT("(0,2) + 1*C^+(1)",      2.0*I + Cd);
 EXPECT_PRINT("(0,2) + 1*A(0)",        2.0*I + A);
 EXPECT_PRINT("(0,2) + sin(t)*A^+(1)", 2.0*I + Ad);

 EXPECT_PRINT("2*t + t^2*C(0)",          "2*t"_te + C);
 EXPECT_PRINT("2*t + 1*C^+(1)",          "2*t"_te + Cd);
 EXPECT_PRINT("2*t + 1*A(0)",            "2*t"_te + A);
 EXPECT_PRINT("2*t + sin(t)*A^+(1)",     "2*t"_te + Ad);
 EXPECT_PRINT("(0,2*t) + t^2*C(0)",      "2*t"_te*I + C);
 EXPECT_PRINT("(0,2*t) + 1*C^+(1)",      "2*t"_te*I + Cd);
 EXPECT_PRINT("(0,2*t) + 1*A(0)",        "2*t"_te*I + A);
 EXPECT_PRINT("(0,2*t) + sin(t)*A^+(1)", "2*t"_te*I + Ad);

 EXPECT_PRINT("1*C^+(1) + t^2*C(0) + sin(t)*A^+(1) + 1*A(0)", C + Cd + A + Ad);

 // Subtraction
 EXPECT_PRINT("-2 + t^2*C(0)",          C - 2.0);
 EXPECT_PRINT("-2 + 1*C^+(1)",          Cd - 2.0);
 EXPECT_PRINT("-2 + 1*A(0)",            A - 2.0);
 EXPECT_PRINT("-2 + sin(t)*A^+(1)",     Ad - 2.0);
 EXPECT_PRINT("(0,-2) + t^2*C(0)",      C - 2.0*I);
 EXPECT_PRINT("(0,-2) + 1*C^+(1)",      Cd - 2.0*I);
 EXPECT_PRINT("(0,-2) + 1*A(0)",        A - 2.0*I);
 EXPECT_PRINT("(0,-2) + sin(t)*A^+(1)", Ad - 2.0*I);

 EXPECT_PRINT("-(2*t) + t^2*C(0)",          C - "2*t"_te);
 EXPECT_PRINT("-(2*t) + 1*C^+(1)",          Cd - "2*t"_te);
 EXPECT_PRINT("-(2*t) + 1*A(0)",            A - "2*t"_te);
 EXPECT_PRINT("-(2*t) + sin(t)*A^+(1)",     Ad - "2*t"_te);
 EXPECT_PRINT("(0,-(2*t)) + t^2*C(0)",      C - "2*t"_te*I);
 EXPECT_PRINT("(0,-(2*t)) + 1*C^+(1)",      Cd - "2*t"_te*I);
 EXPECT_PRINT("(0,-(2*t)) + 1*A(0)",        A - "2*t"_te*I);
 EXPECT_PRINT("(0,-(2*t)) + sin(t)*A^+(1)", Ad - "2*t"_te*I);

 EXPECT_PRINT("2 + -(t^2)*C(0)",          2.0 - C);
 EXPECT_PRINT("2 + -1*C^+(1)",            2.0 - Cd);
 EXPECT_PRINT("2 + -1*A(0)",              2.0 - A);
 EXPECT_PRINT("2 + -(sin(t))*A^+(1)",     2.0 - Ad);
 EXPECT_PRINT("(0,2) + -(t^2)*C(0)",      2.0*I - C);
 EXPECT_PRINT("(0,2) + -1*C^+(1)",        2.0*I - Cd);
 EXPECT_PRINT("(0,2) + -1*A(0)",          2.0*I - A);
 EXPECT_PRINT("(0,2) + -(sin(t))*A^+(1)", 2.0*I - Ad);

 EXPECT_PRINT("2*t + -(t^2)*C(0)",          "2*t"_te - C);
 EXPECT_PRINT("2*t + -1*C^+(1)",            "2*t"_te - Cd);
 EXPECT_PRINT("2*t + -1*A(0)",              "2*t"_te - A);
 EXPECT_PRINT("2*t + -(sin(t))*A^+(1)",     "2*t"_te - Ad);
 EXPECT_PRINT("(0,2*t) + -(t^2)*C(0)",      "2*t"_te*I - C);
 EXPECT_PRINT("(0,2*t) + -1*C^+(1)",        "2*t"_te*I - Cd);
 EXPECT_PRINT("(0,2*t) + -1*A(0)",          "2*t"_te*I - A);
 EXPECT_PRINT("(0,2*t) + -(sin(t))*A^+(1)", "2*t"_te*I - Ad);

 EXPECT_PRINT("-1*C^+(1) + t^2*C(0) + -(sin(t))*A^+(1) + -1*A(0)", C - Cd - A - Ad);

 // Multiplication
 EXPECT_PRINT("(t^2)*(3)*C(0)",          C * 3.0);
 EXPECT_PRINT("3*C^+(1)",                Cd * 3.0);
 EXPECT_PRINT("3*A(0)",                  A * 3.0);
 EXPECT_PRINT("(sin(t))*(3)*A^+(1)",     Ad * 3.0);
 EXPECT_PRINT("(0,(t^2)*(3))*C(0)",      C * 3.0*I);
 EXPECT_PRINT("(0,3)*C^+(1)",            Cd * 3.0*I);
 EXPECT_PRINT("(0,3)*A(0)",              A * 3.0*I);
 EXPECT_PRINT("(0,(sin(t))*(3))*A^+(1)", Ad * 3.0*I);

 EXPECT_PRINT("(t^2)*(2*t)*C(0)",          C * "2*t"_te);
 EXPECT_PRINT("2*t*C^+(1)",                Cd * "2*t"_te);
 EXPECT_PRINT("2*t*A(0)",                  A * "2*t"_te);
 EXPECT_PRINT("(sin(t))*(2*t)*A^+(1)",     Ad * "2*t"_te);
 EXPECT_PRINT("(0,(t^2)*(2*t))*C(0)",      C * "2*t"_te*I);
 EXPECT_PRINT("(0,2*t)*C^+(1)",            Cd * "2*t"_te*I);
 EXPECT_PRINT("(0,2*t)*A(0)",              A * "2*t"_te*I);
 EXPECT_PRINT("(0,(sin(t))*(2*t))*A^+(1)", Ad * "2*t"_te*I);

 EXPECT_PRINT("(t^2)*(3)*C(0)",          3.0 * C);
 EXPECT_PRINT("3*C^+(1)",                3.0 * Cd);
 EXPECT_PRINT("3*A(0)",                  3.0 * A);
 EXPECT_PRINT("(sin(t))*(3)*A^+(1)",     3.0 * Ad);
 EXPECT_PRINT("(0,(t^2)*(3))*C(0)",      3.0*I * C);
 EXPECT_PRINT("(0,3)*C^+(1)",            3.0*I * Cd);
 EXPECT_PRINT("(0,3)*A(0)",              3.0*I * A);
 EXPECT_PRINT("(0,(sin(t))*(3))*A^+(1)", 3.0*I * Ad);

 EXPECT_PRINT("(t^2)*(2*t)*C(0)",          "2*t"_te * C);
 EXPECT_PRINT("2*t*C^+(1)",                "2*t"_te * Cd);
 EXPECT_PRINT("2*t*A(0)",                  "2*t"_te * A);
 EXPECT_PRINT("(sin(t))*(2*t)*A^+(1)",     "2*t"_te * Ad);
 EXPECT_PRINT("(0,(t^2)*(2*t))*C(0)",      "2*t"_te*I * C);
 EXPECT_PRINT("(0,2*t)*C^+(1)",            "2*t"_te*I * Cd);
 EXPECT_PRINT("(0,2*t)*A(0)",              "2*t"_te*I * A);
 EXPECT_PRINT("(0,(sin(t))*(2*t))*A^+(1)", "2*t"_te*I * Ad);

 EXPECT_PRINT("-2 + 2*C^+(2)C(2) + -2*C^+(2)A^+(2) + 2*C(2)A(2) + -1*[A^+(2)]^2 + 1*[A(2)]^2",
              (c<te>(2) + c_dag<te>(2) + a<te>(2) + a_dag<te>(2))*
              (c<te>(2) - c_dag<te>(2) + a<te>(2) - a_dag<te>(2)));

 // (n_up * a + n_dn * a^+)^3
 auto expr = n<te>("up") * a<te>(0) + n<te>("dn") * a_dag<te>(0);
 EXPECT_PRINT("1*C^+(dn)C(dn)A^+(0) + 1*C^+(up)C(up)A(0)", expr);
 expr = expr*expr*expr;
 EXPECT_PRINT("3*C^+(dn)C^+(up)C(up)C(dn)A^+(0) + 3*C^+(dn)C^+(up)C(up)C(dn)A(0) + "
              "1*C^+(dn)C(dn)[A^+(0)]^3 + 1*C^+(up)C(up)[A(0)]^3 + "
              "3*C^+(dn)C^+(up)C(up)C(dn)[A^+(0)]^2A(0) + 3*C^+(dn)C^+(up)C(up)C(dn)A^+(0)[A(0)]^2",
              expr);

 // Dagger
 auto X = te("t^2","sin(t)")*c_dag<te>(1) * c_dag<te>(2) * c<te>(3) * c<te>(4) * a_dag<te>(5) * a<te>(6);
 EXPECT_PRINT("(-(t^2),-(sin(t)))*C^+(1)C^+(2)C(4)C(3)A^+(5)A(6)", X);
 EXPECT_PRINT("(-(t^2),-(-(sin(t))))*C^+(3)C^+(4)C(2)C(1)A^+(6)A(5)", dagger(X));
}

MAKE_MAIN;
