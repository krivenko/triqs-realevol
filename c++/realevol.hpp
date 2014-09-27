#pragma once

#include <set>
#include <string>
#include <utility>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <triqs/parameters.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/h5/map.hpp>
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/draft/hilbert_space_tools/hilbert_space.hpp>
#include <triqs/draft/hilbert_space_tools/state.hpp>

#include "mesh_container.hpp"
#include "time_expr.hpp"
#include "c_time_expr.hpp"

namespace realevol {

using namespace triqs::utility;
using parameters_t = triqs::params::parameters;

using dcomplex = std::complex<double>;
template<typename Value> using dict_t = std::map<std::string,Value>;

template<typename Mesh, bool ComplexOperators = false>
class solver {

    fundamental_operator_set fops;
    hilbert_space hs;
    state<hilbert_space,dcomplex,false> init_state;

    dict_t<mesh_container<double,Mesh>> results;

    boost::mpi::communicator comm;      // define the communicator, here MPI_COMM_WORLD

public:

    using operator_t = typename std::conditional<
        ComplexOperators,
        many_body_operator<c_time_expr>,
        many_body_operator<time_expr>
    >::type;

    solver(std::set<std::string> const& operator_indices);

    void solve(operator_t h, parameters_t params, dict_t<operator_t> observables = {});

    decltype(init_state) & psi0() { return init_state; }

    dcomplex & psi0(std::set<std::string> const& indices) {
        std::set<fundamental_operator_set::indices_t> tmp;
        for(auto const & i : indices) tmp.insert({i});
        return init_state(hs.get_fock_state(fops,tmp));
    }

    decltype(results) const& get_results() const { return results; }
    static parameters_t solve_parameters();
    static void help();

private:

    static void fill_fops(fundamental_operator_set & fops, operator_t const& op);
};

}