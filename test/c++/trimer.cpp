#include <string>
#include <array>
#include <triqs/gfs/meshes/segment.hpp>

#include "../c++/time_expr_r.hpp"
#include "../c++/realevol.hpp"

using namespace realevol;

double hbar = 1.0;
double U = 1.0;
double mu = 0.5*U;
time_expr_r V = "0.3*(1 - exp(-4*t))";

using operator_t = many_body_operator<time_expr_r>;

int main()
{
    std::array<std::string,3> atoms {"1","2","3"};
    std::array<std::string,2> spins {"up","dn"};

    std::set<std::string> all_indices;
    for(auto a : atoms)
        for(auto s : spins)
            all_indices.insert(a+"-"+s);

    auto H = operator_t();

    // Chemical potential
    for(auto a : atoms)
        for(auto s : spins)
            H += -mu*n<time_expr_r>(a+"-"+s);
    // Hubbard interaction
    for(auto a : atoms) H += U*n<time_expr_r>(a+"-up")*n<time_expr_r>(a+"-dn");
    // Hopping between atoms
    for(auto s : spins){
        for(auto a1 : atoms)
        for(auto a2 : atoms){
            if(a1==a2) continue;
            H += V*c_dag<time_expr_r>(a1+"-"+s)*c<time_expr_r>(a2+"-"+s);
        }
    }

    // Time mesh
    triqs::gfs::segment_mesh mesh(0,50.0,501);

    // Solver object
    using solver_t = solver<triqs::gfs::segment_mesh>;
    solver_t S(all_indices);

    // Parameters
    auto params = solver_t::solve_parameters();

    // Observables
    dict_t<operator_t> observables;
    observables["unity"] = operator_t() + 1.0; // ugly

    // psi0: Sz=1/2
    S.psi0({"1-up","2-up","3-dn"}) = 1.0;
    S.solve(H,params,observables);
    // TODO

    // psi0: Sz=0
    S.psi0().amplitudes()() = .0;
    S.psi0({{"1-up"},{"1-dn"},{"2-up"},{"3-dn"}}) = 1.0;
    S.solve(H,params,observables);
    // TODO

    return EXIT_SUCCESS;
}