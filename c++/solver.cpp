#include <triqs/utility/first_include.hpp>

#include <boost/timer/timer.hpp>
#include <boost/variant/apply_visitor.hpp>

#include "solver.hpp"

namespace realevol {

using triqs::hilbert_space::imperative_operator;
using triqs::hilbert_space::sub_hilbert_space;
using triqs::hilbert_space::fock_state_t;

solver::solver(std::set<indices_t> const& operator_indices) :
    fops(operator_indices), hs(fops), init_state(hs) {
}

void solver::solve(solve_parameters_t const& p) {

    // FIXME
    if(comm.size()>1) TRIQS_RUNTIME_ERROR << "Running on more than one MPI node is not yet supported.";

    _last_solve_parameters.reset(new solve_parameters_t(p));

    boost::timer::auto_cpu_timer solve_timer(p.verbosity >= 1 ? "Simulation took %w seconds\n" : "");

    auto sim = simulation(comm,fops,hs,init_state,p);
    results = boost::apply_visitor(sim, p.mesh);
}

}
