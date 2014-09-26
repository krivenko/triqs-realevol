#include <string>
#include <array>
#include <triqs/gfs/meshes/segment.hpp>

#include "../c++/time_expr.hpp"
#include "../c++/realevol.hpp"

using namespace realevol;

double hbar = 1.0;
double U = 1.0;
double mu = 0.5*U;
time_expr V = "0.3*(1 - exp(-4*t))";

using operator_t = many_body_operator<time_expr>;

int main()
{
    std::array<std::string,3> atoms {"1","2","3"};
    std::array<std::string,2> spins {"up","dn"};

    auto H = operator_t();

    // Chemical potential
    for(auto a : atoms)
        for(auto s : spins)
            H += -mu*n<time_expr>(a+"-"+s);
    // Hubbard interaction
    for(auto a : atoms) H += U*n<time_expr>(a+"-up")*n<time_expr>(a+"-dn");
    // Hopping between atoms
    for(auto s : spins){
        for(auto a1 : atoms)
        for(auto a2 : atoms){
            if(a1==a2) continue;
            H += V*c_dag<time_expr>(a1+"-"+s)*c<time_expr>(a2+"-"+s);
        }
    }

    // Time mesh
    triqs::gfs::segment_mesh mesh(0,50.0,501);

    // TODO

    return EXIT_SUCCESS;
}