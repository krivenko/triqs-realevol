/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <type_traits>

#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>

#include "common.hpp"

namespace realevol {

// Hamiltonian interpolation between time slices
enum h_approx {Rectangle, Trapezoid, Simpson} approx;

template<h_approx Approx>
class propagator {

 op_on_subspace_t h;            // Hamiltonian
 dcomplex h_coeff;              // Hamiltonian prefactor in the exponential

 // Preallocated temporary objects
 state_on_subspace_t state_1d;
 matrix<dcomplex> lapack_workspace;

 state_on_subspace_t apply_h(state_on_subspace_t const& st, double t, double dt) const;

public:

 using time_it_t = gf_mesh<retime>::const_iterator;

 propagator(op_on_subspace_t const& h, sub_hilbert_space const& sp,
            double hbar, h_approx approx,
            int lanczos_min_matrix_size);

 void operator()(state_on_subspace_t & st,
                 time_it_t const& t_start, time_it_t const& t_end) const;

private:

 void propagate_1d(state_on_subspace_t & st,
                   time_it_t t_min, time_it_t const& t_max, bool forward) const;
 void propagate_lapack(state_on_subspace_t & st,
                       time_it_t t_min, time_it_t const& t_max, bool forward) const;
 void propagate_lanczos(state_on_subspace_t & st,
                        time_it_t t_min, time_it_t const& t_max, bool forward) const;

 // Pointer to the propagate_* method selected upon construction
 void (propagator::* propagate)(state_on_subspace_t &, time_it_t, time_it_t const&, bool) const;
};

}
