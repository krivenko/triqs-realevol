#include <triqs/utility/first_include.hpp>

#include <boost/timer/timer.hpp>
#include <boost/variant/apply_visitor.hpp>

#include "solver.hpp"

namespace realevol {

using triqs::utility::imperative_operator;
using triqs::utility::sub_hilbert_space;
using triqs::utility::fock_state_t;

template<bool ComplexOperators>
solver<ComplexOperators>::solver(std::set<indices_t> const& operator_indices) :
    fops(operator_indices), hs(fops), init_state(hs)
{
}

template<bool ComplexOperators>
void solver<ComplexOperators>::solve(solve_parameters_t<operator_t> const& p) {

    _last_solve_parameters.reset(new solve_parameters_t<operator_t>(p));

    boost::timer::auto_cpu_timer solve_timer(p.verbosity >= 1 ? "Simulation took %w seconds\n" : "");

    auto sim = simulation<ComplexOperators>(comm,fops,hs,init_state,p);
    results = boost::apply_visitor(sim, p.mesh);
}

// Explicit instantiations
template class solver<false>;
template class solver<true>;

}