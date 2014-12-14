#pragma once

#include "common.hpp"
#include "mesh_container_cyclic.hpp"
#include "runge_kutta.hpp"
#include "lanczos.hpp"

using namespace triqs::utility;

namespace realevol {

template<bool ComplexOp, typename Mesh, ode_solve_method Method> class schroedinger_worker;

// Runge-Kutta version
template<bool ComplexOp, typename Mesh> class schroedinger_worker<ComplexOp,Mesh,method_runge_kutta> {

public:
    using mesh_t = Mesh;
    using solution_t = mesh_container_cyclic<state_on_subspace_t,mesh_t>;

private:
    // RHS part of Schroedinger's equation
    struct rhs_t {
        op_on_subspace_t<ComplexOp> h;
        double hbar;

        rhs_t(op_on_subspace_t<ComplexOp> const& hamiltonian, double hbar) : h(hamiltonian), hbar(hbar) {}
        state_on_subspace_t operator()(state_on_subspace_t const& psi, double t) const { return h(psi,t)/(1_j*hbar); }
    };

    solution_t & solution;
    rhs_t rhs;
    runge_kutta<solution_t,rhs_t> solver;
    typename solution_t::iterator it1, it2;

public:

    schroedinger_worker(op_on_subspace_t<ComplexOp> const& hamiltonian, solution_t & solution, double hbar) :
        rhs(hamiltonian,hbar), solution(solution), solver(rhs),
        it1(std::begin(solution)), it2(it1 + solution.storage_size())
    {}

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
template<bool ComplexOp, typename Mesh> class schroedinger_worker<ComplexOp,Mesh,method_lanczos> {

public:
    using mesh_t = Mesh;
    using solution_t = mesh_container_cyclic<state_on_subspace_t,mesh_t>;

private:
    // RHS part of Schroedinger's equation
    struct rhs_t {
        op_on_subspace_t<ComplexOp> h;

        rhs_t(op_on_subspace_t<ComplexOp> const& hamiltonian) : h(hamiltonian) {}
        state_on_subspace_t operator()(state_on_subspace_t const& psi, double t) const { return h(psi,t); }
    };

    rhs_t rhs;
    solution_t & solution;
    double hbar;
    lanczos<solution_t,rhs_t> solver;
    typename solution_t::iterator it1, it2;

public:

    schroedinger_worker(op_on_subspace_t<ComplexOp> const& hamiltonian, solution_t & solution, double hbar) :
        rhs(hamiltonian), solution(solution), hbar(hbar), solver(rhs),
        it1(std::begin(solution)), it2(it1 + solution.storage_size())
    {}

    // Returns true when the end of the mesh is reached-
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