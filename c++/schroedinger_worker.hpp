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
template<bool ComplexOperators, typename Mesh> class schroedinger_worker<ComplexOperators,Mesh,runge_kutta_method> {
};

// Lanczos version
template<bool ComplexOperators, typename Mesh> class schroedinger_worker<ComplexOperators,Mesh,lanczos_method> {
};

// public:
// 
//     using mesh_t = Mesh;
//     using state_t = state<sub_hilbert_space,std::complex<double>,false>;
//     using solution_t = mesh_container_cyclic<state_t,mesh_t>;
// 
// private:
// 
//     // RHS part of Schroedinger's equation
//     struct rhs_t{
//         imperative_operator<sub_hilbert_space,operator_coeff_t,false> h;
//         double hbar;
//         std::complex<double> I = {0,1};
// 
//         rhs_t(operator_t const& hamiltonian, fundamental_operator_set const& fops, double hbar) :
//             h(hamiltonian, fops), hbar(hbar) {}
// 
//         state_t operator()(state_t const& psi, double t){
//             if(Method==ode_solve_lanczos)
//                 return h(psi,t);
//             else
//                 return h(psi,t)/(I*hbar);
//         }
//     };
// 
//     solution_t solution;
//     rhs_t rhs;
// 
// public:
// 
//     using solver_t = typename std::conditional<
//         Method==ode_solve_runge_kutta,
//         runge_kutta<solution_t,rhs_t>,
//         lanczos<solution_t,rhs_t>
//         >::type;
// 
//     schroedinger_worker(sub_hilbert_space const& hs, fundamental_operator_set const& fops,
//                         operator_t const& hamiltonian, initial_state_t const& psi0,
//                         mesh_t const& mesh, double hbar,
//                         std::size_t stored_psi_values) :
//         solution(mesh,stored_psi_values), rhs(hamiltonian,fops,hbar) {
//     }
// 
// };

}