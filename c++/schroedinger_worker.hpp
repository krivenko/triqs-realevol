#pragma once

#include <complex>
#include <type_traits>

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/draft/hilbert_space_tools/hilbert_space.hpp>
#include <triqs/draft/hilbert_space_tools/state.hpp>
#include <triqs/draft/hilbert_space_tools/imperative_operator.hpp>

#include "solve_parameters.hpp"
#include "time_expr_r.hpp"
#include "time_expr_c.hpp"
#include "mesh_container_cyclic.hpp"
#include "runge_kutta.hpp"
#include "lanczos.hpp"

using namespace triqs::utility;

namespace realevol {

template<bool ComplexOperators, typename Mesh, ode_solve_method Method> class schroedinger_worker;

// Runge-Kutta version
template<bool ComplexOperators, typename Mesh> class schroedinger_worker<ComplexOperators,Mesh,method_runge_kutta> {

public:
    using mesh_t = Mesh;
    using state_t = state<sub_hilbert_space,std::complex<double>,false>;
    using solution_t = mesh_container_cyclic<state_t,mesh_t>;
    using imp_operator_t = imperative_operator<sub_hilbert_space,operator_coeff_t<ComplexOperators>,false>;

private:
    // RHS part of Schroedinger's equation
    struct rhs_t {
        imp_operator_t h;
        double hbar;

        rhs_t(imp_operator_t const& hamiltonian, double hbar) : h(hamiltonian), hbar(hbar) {}
        state_t operator()(state_t const& psi, double t) const { return h(psi,t)/(1_j*hbar); }
    };

    solution_t solution;
    rhs_t rhs;
    runge_kutta<solution_t,rhs_t> solver;
    typename solution_t::iterator it1, it2;

public:

    schroedinger_worker(imp_operator_t const& hamiltonian, state_t const& psi0,
                        mesh_t const& mesh, double hbar,
                        int stored_psi_values) :
        solution(mesh,stored_psi_values), rhs(hamiltonian,hbar), solver(rhs),
        it1(std::begin(solution)), it2(it1 + solution.storage_size())
    {
        solution[0] = psi0;
    }

    // Returns true when the end of the mesh is reached
    bool operator()() {
        solver(it1,it2);
        if(it2 == std::end(solution)) return true;
        it1 += solution.storage_size()-1;
        it2 += solution.storage_size()-1;
        if(it2 > std::end(solution)) it2 = std::end(solution);
        return false;
    }
};

// Lanczos version
template<bool ComplexOperators, typename Mesh> class schroedinger_worker<ComplexOperators,Mesh,method_lanczos> {

public:
    using mesh_t = Mesh;
    using state_t = state<sub_hilbert_space,std::complex<double>,false>;
    using solution_t = mesh_container_cyclic<state_t,mesh_t>;
    using imp_operator_t = imperative_operator<sub_hilbert_space,operator_coeff_t<ComplexOperators>,false>;

private:
    // RHS part of Schroedinger's equation
    struct rhs_t {
        imp_operator_t h;

        rhs_t(imp_operator_t const& hamiltonian) : h(hamiltonian) {}
        state_t operator()(state_t const& psi, double t) const { return h(psi,t); }
    };

    solution_t solution;
    rhs_t rhs;
    double hbar;
    lanczos<solution_t,rhs_t> solver;
    typename solution_t::iterator it1, it2;

public:

    schroedinger_worker(imp_operator_t const& hamiltonian, state_t const& psi0,
                        mesh_t const& mesh, double hbar,
                        int stored_psi_values) :
        solution(mesh,stored_psi_values), rhs(hamiltonian), hbar(hbar), solver(rhs),
        it1(std::begin(solution)), it2(it1 + solution.storage_size())
    {
        solution[0] = psi0;
    }

    // Returns true when the end of the mesh is reached
    bool operator()() {
        solver(it1,it2,1.0/(1_j*hbar));
        if(it2 == std::end(solution)) return true;
        it1 += solution.storage_size()-1;
        it2 += solution.storage_size()-1;
        if(it2 > std::end(solution)) it2 = std::end(solution);
        return false;
    }
};

}