#pragma once

#include <string>
#include <utility>
#include <boost/mpi/communicator.hpp>
#include <triqs/parameters.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/h5/map.hpp>

#include "mesh_container.hpp"
#include "time_expr.hpp"
#include "callable_complex.hpp"

namespace realevol {

using triqs::utility::many_body_operator;
using parameters_t = triqs::params::parameters;
template<typename Value> using dict_t = std::map<std::string,Value>;

template<typename Mesh, bool ComplexOperators = false>
class solver {

    using results_t = dict_t<mesh_container<Mesh,double>>;

    results_t _results;
    boost::mpi::communicator comm;      // define the communicator, here MPI_COMM_WORLD

public:

    using operator_t = std::conditional<
        ComplexOperators,
        many_body_operator<callable_complex<time_expr>>,
        many_body_operator<time_expr>
    >;

    solver();

    void solve(operator_t h, parameters_t params, dict_t<operator_t> observables = {});

    results_t const& results() const { return _results; }

    static parameters_t solve_parameters();
    static void help();
};

}