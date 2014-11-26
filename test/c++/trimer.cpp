#include <triqs/utility/first_include.hpp>

#include <string>
#include <array>
#include <triqs/gfs/meshes/segment.hpp>

#include "../c++/time_expr_r.hpp"
#include "../c++/solver.hpp"

using namespace realevol;

using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;

double hbar = 1.0;
double U = 1.0;
double mu = 0.5*U;
time_expr_r V = "0.3*(1 - exp(-4*t))";

// Real-valued version of the solver
using solver_t = solver<false>;

int main()
{
    std::array<int,3> atoms {1,2,3};
    std::array<std::string,2> spins {"up","dn"};

    std::set<solver_t::indices_t> all_indices;
    for(auto a : atoms)
        for(auto s : spins)
            all_indices.insert({a,s});

    auto H = operator_t<false>();

    // Chemical potential
    for(auto a : atoms)
        for(auto s : spins)
            H += -mu*n<time_expr_r>(a,s);
    // Hubbard interaction
    for(auto a : atoms) H += U*n<time_expr_r>(a,"up")*n<time_expr_r>(a,"dn");
    // Hopping between atoms
    for(auto s : spins){
        for(auto a1 : atoms)
        for(auto a2 : atoms){
            if(a1==a2) continue;
            H += V*c_dag<time_expr_r>(a1,s)*c<time_expr_r>(a2,s);
        }
    }

    // Time mesh
    triqs::gfs::segment_mesh mesh(0,50.0,501);

    // Solver object
    solver_t S(all_indices);

    // Parameters
    auto params = solve_parameters_t<false>(H,mesh);
    params.verbosity = 2;
    params.stored_psi_values = 20;

    // Observables
    std::map<std::string,operator_t<false>> observables;
    params.observables["unity"] = operator_t<false>() + 1.0; // ugly

    // psi0: Sz=1/2
    S.psi0({{1,"up"},{2,"up"},{3,"dn"}}) = 1.0;
    params.method = method_lanczos;
    S.solve(params);
    // TODO

    // psi0: Sz=0
    S.psi0().amplitudes()() = .0;
    S.psi0({{1,"up"},{1,"dn"},{2,"up"},{3,"dn"}}) = 1.0;
    params.method = method_runge_kutta;
    S.solve(params);
    // TODO

    return EXIT_SUCCESS;
}