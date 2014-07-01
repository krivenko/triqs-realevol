#include <cstdlib>
#include <cmath>
#include <complex>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "time_expr.hpp"
#include "uniform_mesh.hpp"

using namespace realevol;

int main()
{   
    {
    // Test real-valued expressions
    using time_expr = time_expr<false>;
    #include "time_expr.body.hpp"
    }
    {
    // Test complex-valued expressions
    using time_expr = time_expr<true>;
    #include "time_expr.body.hpp"
    }
    return EXIT_SUCCESS;
}