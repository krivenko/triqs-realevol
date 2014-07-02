#include <triqs/utility/first_include.hpp>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "many_body_operator.hpp"
#include "time_expr.hpp"

using namespace realevol;

template<typename op>
void test_commutators(std::vector<op> const& Cd, std::vector<op> const& C)
{
    std::cout << std::endl << "Commutators:" << std::endl;
    for(auto const& cdi : Cd)
    for(auto const& ci : C){
        std::cout << "[" << cdi << ", " << ci << "] = " << cdi*ci - ci*cdi << std::endl;
    }    
}

template<typename op>
void test_anticommutators(std::vector<op> const& Cd, std::vector<op> const& C)
{
    std::cout << std::endl << "Anticommutators:" << std::endl;
    for(auto const& cdi : Cd)
    for(auto const& ci : C){
        std::cout << "{" << cdi << ", " << ci << "} = " << cdi*ci + ci*cdi << std::endl;
    }
}

template<typename op>
void test_serialization(op const& X, std::string const& name)
{
    // Serialization
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa & X;

    boost::archive::text_iarchive ia(ss);
    op new_X;
    ia & new_X;
    
    std::cout << "New " << name << "= " << new_X << std::endl;
}

int main()
{
    
    {
    // Test real-valued operators
    using scalar_t = double;
    using operator_t = many_body_operator<scalar_t>;
    
    std::cout << std::endl;
    std::cout << "I. Real-valued operators" << std::endl;
    std::cout << "========================" << std::endl;
    
    // Operators without indices
    operator_t op_with_no_indices = c<scalar_t>() + c_dag<scalar_t>() - n<scalar_t>();
    std::cout << "op_with_no_indices = " << op_with_no_indices << std::endl;
    
    // Operators with many indices
    auto op_with_many_indices = c<scalar_t>(1,0.2,"a",true,-2) +
                                c_dag<scalar_t>(3,0.15,"b",false,-5);
    std::cout << "op_with_many_indices = " << op_with_many_indices << std::endl;
                            
    // Commutation relations
    std::vector<operator_t> C = {c<scalar_t>(1), c<scalar_t>(2), c<scalar_t>(3)};
    std::vector<operator_t> Cd = {c_dag<scalar_t>(1), c_dag<scalar_t>(2), c_dag<scalar_t>(3)};
    test_anticommutators(Cd,C);
    test_commutators(Cd,C);
    
    // Algebra
    auto x = c<scalar_t>(0);
    auto y = c_dag<scalar_t>(1);

    std::cout << std::endl << "Algebra:" << std::endl;    
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;

    std::cout << "-x = " << -x << std::endl;
    std::cout << "x + 2.0 = " << x + 2.0 << std::endl;
    std::cout << "2.0 + x = " << 2.0 + x << std::endl;
    std::cout << "x - 2.0 = " << x - 2.0 << std::endl;
    std::cout << "2.0 - x = " << 2.0 - x << std::endl;
    std::cout << "3.0*y = " << 3.0*y << std::endl;
    std::cout << "y*3.0 = " << y*3.0 << std::endl;
    std::cout << "x + y = " << x + y << std::endl;
    std::cout << "x - y = " << x - y << std::endl;
    std::cout << "(x + y)*(x - y) = " << (x + y)*(x - y) << std::endl;

    // N^3
    std::cout << std::endl << "N^3:" << std::endl;
    auto N = n<scalar_t>("up") + n<scalar_t>("dn");
    auto N3 = N*N*N;
    std::cout << "N = " << N << std::endl;
    std::cout << "N^3 = " << N3 << std::endl;
    
    // Serialization
    test_serialization(N3,"N^3");
   
    auto X = c_dag<scalar_t>(1) * c_dag<scalar_t>(2) * c<scalar_t>(3) * c<scalar_t>(4);
    std::cout  << "X = "<< X<<std::endl; 
    std::cout  << "dagger(X) = "<< dagger(X)<<std::endl; 
    }
    
    {
    // Test complex-valued operators
    using scalar_t = std::complex<double>;
    using operator_t = many_body_operator<scalar_t>;
    
    std::cout << std::endl;    
    std::cout << "II. Complex-valued operators" << std::endl;
    std::cout << "============================" << std::endl;
    
    // Operators without indices
    operator_t op_with_no_indices = c<scalar_t>() + c_dag<scalar_t>() - n<scalar_t>();
    std::cout << "op_with_no_indices = " << op_with_no_indices << std::endl;
    
    // Operators with many indices
    auto op_with_many_indices = c<scalar_t>(1,0.2,"a",true,-2) +
                                c_dag<scalar_t>(3,0.15,"b",false,-5);
    std::cout << "op_with_many_indices = " << op_with_many_indices << std::endl;
    
    // Commutation relations
    std::vector<operator_t> C = {c<scalar_t>(1), c<scalar_t>(2), c<scalar_t>(3)};
    std::vector<operator_t> Cd = {c_dag<scalar_t>(1), c_dag<scalar_t>(2), c_dag<scalar_t>(3)};
    test_anticommutators(Cd,C);
    test_commutators(Cd,C);
    
    // Algebra
    auto x = c<scalar_t>(0);
    auto y = c_dag<scalar_t>(1);

    auto I = scalar_t{0,1.0};
    
    std::cout << std::endl << "Algebra:" << std::endl;    
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;

    std::cout << "-x = " << -x << std::endl;
    std::cout << "x + 2.0*I = " << x + 2.0*I << std::endl;
    std::cout << "2.0 + x = " << 2.0 + x << std::endl;
    std::cout << "x - 2.0 = " << x - 2.0 << std::endl;
    std::cout << "2.0*I - x = " << 2.0*I - x << std::endl;
    std::cout << "3.0*y = " << 3.0*y << std::endl;
    std::cout << "y*3.0*I = " << y*3.0*I << std::endl;
    std::cout << "x + y = " << x + y << std::endl;
    std::cout << "x - y = " << x - y << std::endl;
    std::cout << "(x + y)*(x - y) = " << (x + y)*(x - y) << std::endl;

    // N^3
    std::cout << std::endl << "N^3:" << std::endl;
    auto N = n<scalar_t>("up") + n<scalar_t>("dn");
    auto N3 = N*N*N;
    std::cout << "N = " << N << std::endl;
    std::cout << "N^3 = " << N3 << std::endl;
    
    // Serialization
    test_serialization(N3,"N^3");
   
    auto X = I*c_dag<scalar_t>(1) * c_dag<scalar_t>(2) * c<scalar_t>(3) * c<scalar_t>(4);
    std::cout  << "X = "<< X<<std::endl; 
    std::cout  << "dagger(X) = "<< dagger(X)<<std::endl; 
    
    }
    
    /*{
    // Test real-expression-valued operators
    using scalar_t = time_expr<false>;
    using operator_t = many_body_operator<scalar_t>;

    std::cout << std::endl;
    std::cout << "III. Real-expression-valued operators" << std::endl;
    std::cout << "=====================================" << std::endl;
    
    // Operators without indices
    operator_t op_with_no_indices = c<scalar_t>() + c_dag<scalar_t>() - n<scalar_t>();
    std::cout << "op_with_no_indices = " << op_with_no_indices << std::endl;
    
    // Operators with many indices
    auto op_with_many_indices = c<scalar_t>(1,0.2,"a",true,-2) +
                                c_dag<scalar_t>(3,0.15,"b",false,-5);
    std::cout << "op_with_many_indices = " << op_with_many_indices << std::endl;
    }*/
    
    return EXIT_SUCCESS;
}