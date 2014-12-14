#pragma once

#include <set>
#include <string>
#include <utility>
#include <memory>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <triqs/h5/map.hpp>

#include "common.hpp"
#include "any_mesh.hpp"
#include "mesh_container.hpp"
#include "solve_parameters.hpp"
#include "simulation.hpp"

namespace realevol {

using dcomplex = std::complex<double>;

template<bool ComplexOp = false>
class solver {

    fundamental_operator_set fops;
    hilbert_space hs;
    state<hilbert_space,dcomplex,false> init_state;

    var_results_t results;

    boost::mpi::environment env;
    boost::mpi::communicator comm;      // define the communicator, here MPI_COMM_WORLD

public:

    using indices_t = typename operator_t<ComplexOp>::indices_t;

    solver(std::set<indices_t> const& operator_indices);

    TRIQS_WRAP_ARG_AS_DICT
    void solve(solve_parameters_t<ComplexOp> const& p);

    /// Set of parameters used in the last call to solve
    solve_parameters_t<ComplexOp> get_last_solve_parameters() const {return *_last_solve_parameters;}

    decltype(init_state) & psi0() { return init_state; }
    dcomplex & psi0(std::set<indices_t> const& indices) { return init_state(hs.get_fock_state(fops,indices)); }

    decltype(results) const& get_results() const { return results; }

private:

    std::unique_ptr<solve_parameters_t<ComplexOp>> _last_solve_parameters; // parameters of the last call to solve
};

}