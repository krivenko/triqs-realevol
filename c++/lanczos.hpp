#pragma once

#include <cmath>
#include "lanczos_worker.hpp"

namespace realevol {

template<
    class SolutionMeshContainer,
    class RHSFunction
> struct lanczos {

    using solution_mesh_container_t = SolutionMeshContainer;
    using mesh_t = typename solution_mesh_container_t::mesh_t;
    using var_t = typename mesh_t::domain_pt_t;
    using value_t = typename solution_mesh_container_t::value_type;
    using scalar_t = typename value_t::value_type;
    using rhs_t = RHSFunction;

    // The right-hand part must be representable in a form
    // rhs_prefactor * rhs, where rhs is a Hermitian operator.
    lanczos(const rhs_t& hermitian_rhs /* must be Hermitian */, double gs_energy_convergence = 1e-10) :
        hermitian_rhs(hermitian_rhs), gs_energy_convergence(gs_energy_convergence) {}

    template<typename Iterator>
    void operator()(Iterator first, Iterator last, scalar_t rhs_prefactor = 1.0)
    // Takes the initial value from *first
    // Fills a range (first,last)
    {
        auto next = first; ++next;

        // Value of the independent variable
        var_t x;

        // RHS as a function of the solution only
        using std::placeholders::_1;
        auto rhs_of_sol_only = std::bind(hermitian_rhs, _1, std::cref(x));

        // Lanczos worker object
        lanczos_worker<decltype(rhs_of_sol_only),value_t> lw(rhs_of_sol_only,gs_energy_convergence);

        // Dimension of the problem
        std::size_t N = (first->value).size();

        triqs::arrays::matrix<scalar_t> lanczos_exp(N,N);

        for(;next != last; ++first, ++next){
            // Value of the independent variable
            x = 0.5*(first->mesh_point + next->mesh_point);

            // Solution at current point
            value_t const& U(first->value);

            // Construct the Krylov basis
            auto norm = std::sqrt(dot_product(U,U));
            lw(U/norm);

            // Step of the mesh
            var_t step = next->mesh_point - first->mesh_point;

            // Construct the propagation exponent
            auto eigenvalues = lw.values();
            std::size_t krylov_dim = eigenvalues.size();
            auto all = range(0,krylov_dim);

            for (std::size_t n = 0; n < krylov_dim; ++n)
                lanczos_exp(n,all) = exp(step * rhs_prefactor * eigenvalues(n)) * lw.vectors()(n,all);
            lanczos_exp(all,all) = lw.vectors().transpose() * lanczos_exp(all,all);

            // Propagate
            auto krylov_coeffs = norm * lanczos_exp(all, 0);
            next->value = lw.krylov_2_fock(krylov_coeffs);
        }
    }

private:
    const rhs_t& hermitian_rhs;
    const double gs_energy_convergence;
};

}