#pragma once

#include<functional>

namespace realevol {

template<
    class SolutionMeshContainer,
    class RHSFunction
> struct runge_kutta {

    using solution_mesh_container_t = SolutionMeshContainer;
    using mesh_t = typename solution_mesh_container_t::mesh_t;
    using var_t = typename mesh_t::domain_pt_t;
    using value_t = typename solution_mesh_container_t::value_type;
    using rhs_t = RHSFunction;

    runge_kutta(const rhs_t& rhs) : rhs(rhs) {}

    template<typename Iterator>
    void operator()(Iterator first, Iterator last)
    // Takes the initial value from *first
    // Fills a range (first,last)
    {
        auto next = first; ++next;
        // Step of a uniform mesh
        var_t step = next->mesh_point - first->mesh_point;

        // Number of unknown functions
        std::size_t N = (first->value).size();

        // Temporaries
        value_t rhs_arg(N), X1(N), X2(N), X3(N), X4(N);

        using std::placeholders::_1;
        auto rhs_of_var_only = std::bind(rhs, std::cref(rhs_arg), _1);

        for(;next != last; ++first, ++next){

            // Value of the independent variable
            var_t x = first->mesh_point;

            // Solution at this point
            value_t const& U(first->value);

            // Stage 1
            rhs_arg = U;
            X1 = rhs_of_var_only(x);

            // Stage 2
            rhs_arg = U + 0.5*step*X1;
            X2 = rhs_of_var_only(x + 0.5*step);

            // Stage 3
            rhs_arg = U + 0.5*step*X2;
            X3 = rhs_of_var_only(x + 0.5*step);

            // Stage 4
            rhs_arg = U + step*X3;
            X4 = rhs_of_var_only(x + step);

            next->value = U + step/6.0*(X1 + 2.0*X2 + 2.0*X3 + X4);
        }
    }

private:
    const rhs_t& rhs;
};

}