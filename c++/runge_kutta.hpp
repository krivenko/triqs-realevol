#pragma once

#include<functional>

namespace realevol {

template<
    class SolutionMeshContainer,
    class RHSFunction
> struct runge_kutta {

    using solution_mesh_container_t = SolutionMeshContainer;
    using mesh_t = typename solution_mesh_container_t::mesh_t;
    using var_t = typename mesh_t::mesh_point_t;
    using value_t = typename solution_mesh_container_t::value_type;
    using rhs_t = RHSFunction;

    runge_kutta(const rhs_t& rhs) : rhs(rhs) {}

    void operator()(solution_mesh_container_t & solution, value_t const& initial_value)
    {
        // Step of a uniform mesh
        var_t step = solution.get_mesh().get_step();

        // Number of unknown functions
        std::size_t N = initial_value.size();

        // Set the initial value
        auto it = std::begin(solution);
        it->value = initial_value;

        // Temporaries
        value_t rhs_arg(N), X1(N), X2(N), X3(N), X4(N);

        using std::placeholders::_1;
        auto rhs_of_var_only = std::bind(rhs, std::cref(rhs_arg), _1);

        while(true){
            // Value of the independent variable
            var_t x = it->mesh_point.value;

            // Solution at this point
            value_t const& U(it->value);

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

            ++it;
            if(it == std::end(solution)) break;

            it->value = U + step/6.0*(X1 + 2.0*X2 + 2.0*X3 + X4);
        }
    }

private:
    const rhs_t& rhs;
};

}