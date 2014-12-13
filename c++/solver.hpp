#pragma once

#include <set>
#include <string>
#include <utility>
#include <memory>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "time_expr_r.hpp"
#include "time_expr_c.hpp"

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/h5/map.hpp>

#include "any_mesh.hpp"
#include "mesh_container.hpp"
#include "time_expr_r.hpp"
#include "time_expr_c.hpp"
#include "solve_parameters.hpp"
#include "simulation.hpp"

namespace realevol {

using dcomplex = std::complex<double>;

template<bool ComplexOperators = false>
class solver {

    fundamental_operator_set fops;
    hilbert_space hs;
    state<hilbert_space,dcomplex,false> init_state;

    results_t results;

    boost::mpi::environment env;
    boost::mpi::communicator comm;      // define the communicator, here MPI_COMM_WORLD

public:

    using indices_t = typename operator_t<ComplexOperators>::indices_t;

    solver(std::set<indices_t> const& operator_indices);

    TRIQS_WRAP_ARG_AS_DICT
    void solve(solve_parameters_t<ComplexOperators> const& p);

    /// Set of parameters used in the last call to solve
    solve_parameters_t<ComplexOperators> get_last_solve_parameters() const {return *_last_solve_parameters;}

    decltype(init_state) & psi0() { return init_state; }
    dcomplex & psi0(std::set<indices_t> const& indices) { return init_state(hs.get_fock_state(fops,indices)); }

    decltype(results) const& get_results() const { return results; }

private:

    std::unique_ptr<solve_parameters_t<ComplexOperators>> _last_solve_parameters; // parameters of the last call to solve
};

}