#pragma once

#include <string>
#include <utility>
#include <boost/mpi/communicator.hpp>
#include <triqs/parameters.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/h5/map.hpp>
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>

#include "mesh_container.hpp"
#include "time_expr.hpp"
#include "c_time_expr.hpp"

namespace realevol {

using namespace triqs::utility;
using parameters_t = triqs::params::parameters;
template<typename Value> using dict_t = std::map<std::string,Value>;

template<typename Mesh, bool ComplexOperators = false>
class solver {

    using results_t = dict_t<mesh_container<Mesh,double>>;

    results_t _results;
    boost::mpi::communicator comm;      // define the communicator, here MPI_COMM_WORLD

public:

    using operator_t = typename std::conditional<
        ComplexOperators,
        many_body_operator<c_time_expr>,
        many_body_operator<time_expr>
    >::type;

    solver();

    void solve(operator_t h, parameters_t params, dict_t<operator_t> observables = {});

    results_t const& results() const { return _results; }

    static parameters_t solve_parameters();
    static void help();

private:

    static void fill_fops(fundamental_operator_set & fops, operator_t const& op);
};

}