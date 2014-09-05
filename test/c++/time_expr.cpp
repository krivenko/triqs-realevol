#include <cstdlib>
#include <cmath>
#include <complex>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "time_expr.hpp"
#include "callable_complex.hpp"
#include "uniform_mesh.hpp"

using namespace realevol;

int main()
{
    {
    // Test real-valued expressions
    using expr_t = time_expr;

    expr_t te1("t^2"), te2("t + sin(_pi/2)"), te3("sqrt(9.0) + 1.5");

    // Check whether the expressions really depend on time
    if(is_constant(te1) || is_constant(te2) || !is_constant(te3))
        return EXIT_FAILURE;

    // Write to a archive
    std::stringstream archive_str;
    boost::archive::text_oarchive oa(archive_str);
    oa << te1; oa << te2; oa << te3;

    // Read from an archive
    boost::archive::text_iarchive ia(archive_str);
    expr_t read_expr;
    ia >> read_expr; if(read_expr != te1) return EXIT_FAILURE;
    ia >> read_expr; if(read_expr != te2) return EXIT_FAILURE;
    ia >> read_expr; if(read_expr != te3) return EXIT_FAILURE;

    // Check correctness of numerical expressions
    double T[] = {0, 0.1, 10, 55};
    double TE1_res[] = {0, 0.01, 100, 3025};
    double TE2_res[] = {1, 1.1, 11, 56};
    double TE3_res[] = {4.5, 4.5, 4.5, 4.5};

    for(int i = 0; i < 4; ++i){
        if(std::abs(te1(T[i]) - TE1_res[i]) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te2(T[i]) - TE2_res[i]) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te3(T[i]) - TE3_res[i]) >= 1e-10) return EXIT_FAILURE;
    }

    // Unary minus
    expr_t mte2 = -te2;
    for(int i = 0; i < 4; ++i){
        if(std::abs(mte2(T[i]) - (-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Addition of expressions
    expr_t te1pte2 = te1 + te2;
    expr_t te1phalf = te1 + 0.5;
    expr_t halfpte2 = 0.5 + te2;

    for(int i = 0; i < 4; ++i){
        if(std::abs(te1pte2(T[i]) - (TE1_res[i]+TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1phalf(T[i]) - (TE1_res[i]+0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfpte2(T[i]) - (0.5+TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Subtraction of expressions
    expr_t te1mte2 = te1 - te2;
    expr_t te1mhalf = te1 - 0.5;
    expr_t halfmte2 = 0.5 - te2;
    for(int i = 0; i < 4; ++i){
        if(std::abs(te1mte2(T[i]) - (TE1_res[i]-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1mhalf(T[i]) - (TE1_res[i]-0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfmte2(T[i]) - (0.5-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Multiplication of expressions
    expr_t te1ppte2 = te1 * te2;
    expr_t te1pphalf = te1 * 0.5;
    expr_t halfppte2 = 0.5 * te2;
    for(int i = 0; i < 4; ++i){
        if(std::abs(te1ppte2(T[i]) - (TE1_res[i]*TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1pphalf(T[i]) - (TE1_res[i]*0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfppte2(T[i]) - (0.5*TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Division of expressions 
    expr_t te1dte2 = te1 / te2;
    expr_t te1dhalf = te1 / 0.5;
    expr_t halfdte2 = 0.5 / te2;
    for(int i = 0; i < 4; ++i){
        if(std::abs(te1dte2(T[i]) - (TE1_res[i]/TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1dhalf(T[i]) - (TE1_res[i]/0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfdte2(T[i]) - (0.5/TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Test try_reduce_to_constant()
    expr_t te0t("0*t"), te1t("1*t");
    if(is_constant(te0t) || is_constant(te1t)) return EXIT_FAILURE;
    uniform_mesh<> m(0,100,101);
    try_reduce_to_constant(te0t, m);
    try_reduce_to_constant(te1t, m);
    if(!is_constant(te0t)) return EXIT_FAILURE;
    if(is_constant(te1t)) return EXIT_FAILURE;

    // Test assignments
    expr_t tea;
    tea = te1t; if(tea != te1t) return EXIT_FAILURE;
    tea = "1*t"; if(tea != te1t) return EXIT_FAILURE;
    tea = std::string("1*t"); if(tea != te1t) return EXIT_FAILURE;
    tea = .0; if(tea != te0t) return EXIT_FAILURE;

    // Test user-defined literals
    expr_t tel1 = "t^3-2*t";
    expr_t tel2 = 0.3;
    if(tel1 != "t^3-2*t"_te) return EXIT_FAILURE;
    if(tel2 != 0.3_te) return EXIT_FAILURE;
    }
    {
    // Test complex-valued expressions
    using expr_t = callable_complex<time_expr>;

    const std::complex<double> I(0,1);

    expr_t te1("t^2",1.0), te2("t + sin(_pi/2)","t^3");
    expr_t te3("sqrt(9.0) + 1.5","2.9"), te4(1.9+I*2.8);

    // Check whether the expressions really depend on time
    if(is_constant(te1) || is_constant(te2) || !is_constant(te3) || !is_constant(te4))
        return EXIT_FAILURE;

    // Write to a archive
    std::stringstream archive_str;
    boost::archive::text_oarchive oa(archive_str);
    oa << te1; oa << te2; oa << te3;

    // Read from an archive
    boost::archive::text_iarchive ia(archive_str);
    expr_t read_expr;
    ia >> read_expr; if(read_expr != te1) return EXIT_FAILURE;
    ia >> read_expr; if(read_expr != te2) return EXIT_FAILURE;
    ia >> read_expr; if(read_expr != te3) return EXIT_FAILURE;

    // Test assignments
    expr_t tea;
    tea = expr_t("t^2",1.0); if(tea != expr_t("t^2",1.0)) return EXIT_FAILURE;
    tea = 0.7; if(tea != expr_t("0.7")) return EXIT_FAILURE;
    tea = std::string("2*t+4*t^3"); if(tea != expr_t("2*t+4*t^3")) return EXIT_FAILURE;
    tea = "8-t"; if(tea != expr_t("8-t")) return EXIT_FAILURE;
    tea = {2.0,5.0}; if(tea != expr_t("2.0","5.0")) return EXIT_FAILURE;
    if(tea != expr_t(2.0,5.0)) return EXIT_FAILURE;

    // Check correctness of numerical expressions
    double T[] = {0, 0.1, 10, 55};
    std::complex<double> TE1_res[] = {{0,1.0}, {0.01,1.0}, {100,1.0}, {3025,1.0}};
    std::complex<double> TE2_res[] = {{1,0}, {1.1,0.001}, {11,1000}, {56,166375}};
    std::complex<double> TE3_res[] = {{4.5,2.9}, {4.5,2.9}, {4.5,2.9}, {4.5,2.9}};
    std::complex<double> TE4_res[] = {{1.9,2.8}, {1.9,2.8}, {1.9,2.8}, {1.9,2.8}};

    for(int i = 0; i < 4; ++i){
        if(std::abs(te1(T[i]) - TE1_res[i]) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te2(T[i]) - TE2_res[i]) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te3(T[i]) - TE3_res[i]) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te4(T[i]) - TE4_res[i]) >= 1e-10) return EXIT_FAILURE;
    }

    // Unary minus
    expr_t mte2 = -te2;
    for(int i = 0; i < 4; ++i){
        if(std::abs(mte2(T[i]) - (-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Addition of expressions
    expr_t te1pte2 = te1 + te2;
    expr_t te1phalf = te1 + 0.5;
    expr_t halfpte2 = 0.5 + te2;
    expr_t te1pihalf = te1 + 0.5*I;
    expr_t ihalfpte2 = 0.5*I + te2;

    for(int i = 0; i < 4; ++i){
        if(std::abs(te1pte2(T[i]) - (TE1_res[i]+TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1phalf(T[i]) - (TE1_res[i]+0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfpte2(T[i]) - (0.5+TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1pihalf(T[i]) - (TE1_res[i]+0.5*I)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(ihalfpte2(T[i]) - (0.5*I+TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Subtraction of expressions
    expr_t te1mte2 = te1 - te2;
    expr_t te1mhalf = te1 - 0.5;
    expr_t halfmte2 = 0.5 - te2;
    expr_t te1mihalf = te1 - 0.5*I;
    expr_t ihalfmte2 = 0.5*I - te2;
    for(int i = 0; i < 4; ++i){
        if(std::abs(te1mte2(T[i]) - (TE1_res[i]-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1mhalf(T[i]) - (TE1_res[i]-0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfmte2(T[i]) - (0.5-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1mihalf(T[i]) - (TE1_res[i]-0.5*I)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(ihalfmte2(T[i]) - (0.5*I-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Multiplication of expressions
    expr_t te1ppte2 = te1 * te2;
    expr_t te1pphalf = te1 * 0.5;
    expr_t halfppte2 = 0.5 * te2;
    expr_t te1ppihalf = te1 * 0.5*I;
    expr_t ihalfppte2 = 0.5*I * te2;
    for(int i = 0; i < 4; ++i){
        if(std::abs(te1ppte2(T[i]) - (TE1_res[i]*TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1pphalf(T[i]) - (TE1_res[i]*0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfppte2(T[i]) - (0.5*TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1ppihalf(T[i]) - (TE1_res[i]*0.5*I)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(ihalfppte2(T[i]) - (0.5*I*TE2_res[i])) >= 1e-10) return EXIT_FAILURE;        
    }

    // Division of expressions 
    expr_t te1dte2 = te1 / te2;
    expr_t te1dhalf = te1 / 0.5;
    expr_t halfdte2 = 0.5 / te2;
    expr_t te1dihalf = te1 / (0.5*I);
    expr_t ihalfdte2 = (0.5*I) / te2;    
    for(int i = 0; i < 4; ++i){
        if(std::abs(te1dte2(T[i]) - (TE1_res[i]/TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1dhalf(T[i]) - (TE1_res[i]/0.5)) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(halfdte2(T[i]) - (0.5/TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(te1dihalf(T[i]) - (TE1_res[i]/(0.5*I))) >= 1e-10) return EXIT_FAILURE;
        if(std::abs(ihalfdte2(T[i]) - ((0.5*I)/TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    }

    // Test try_reduce_to_constant()
    expr_t te0t("0*t","2"), te1t("1*t","t^3");
    if(is_constant(te0t) || is_constant(te1t)) return EXIT_FAILURE;
    uniform_mesh<> m(0,100,101);
    try_reduce_to_constant(te0t, m);
    try_reduce_to_constant(te1t, m);
    if(!is_constant(te0t)) return EXIT_FAILURE;
    if(is_constant(te1t)) return EXIT_FAILURE;

    }
    return EXIT_SUCCESS;
}