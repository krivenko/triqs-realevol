#pragma once

#include <set>
#include <string>
#include <utility>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/h5/map.hpp>

#include "mesh_container.hpp"
#include "time_expr_r.hpp"
#include "time_expr_c.hpp"
#include "solve_parameters.hpp"

#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/draft/hilbert_space_tools/hilbert_space.hpp>
#include <triqs/draft/hilbert_space_tools/state.hpp>

namespace realevol {

using triqs::utility::many_body_operator;
using triqs::utility::fundamental_operator_set;
using triqs::utility::hilbert_space;
using triqs::utility::state;

using dcomplex = std::complex<double>;

template<typename Mesh, bool ComplexOperators = false>
class solver {

    fundamental_operator_set fops;
    hilbert_space hs;
    state<hilbert_space,dcomplex,false> init_state;

    std::map<std::string,mesh_container<double,Mesh>> results;

    boost::mpi::communicator comm;      // define the communicator, here MPI_COMM_WORLD

public:

    using operator_coeff_t = typename std::conditional<ComplexOperators,time_expr_c,time_expr_r>::type;
    using operator_t = many_body_operator<operator_coeff_t>;

    solver(std::set<std::string> const& operator_indices);

    TRIQS_WRAP_ARG_AS_DICT
    void solve(solve_parameters_t<operator_t,Mesh> const& p);

    decltype(init_state) & psi0() { return init_state; }

    dcomplex & psi0(std::set<std::string> const& indices) {
        std::set<fundamental_operator_set::indices_t> tmp;
        for(auto const & i : indices) tmp.insert({i});
        return init_state(hs.get_fock_state(fops,tmp));
    }

    decltype(results) const& get_results() const { return results; }
    static void help();
};

}