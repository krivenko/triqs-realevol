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
#include <memory>

#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>

#include "common.hpp"

namespace realevol {

// Hamiltonian interpolation between time slices
enum h_interpolation {Rectangle, Trapezoid, Simpson};

class propagator;

// 'Diagonalization' on 1-d subspaces
struct diag_1d_t {
 propagator const& prop;
 state_on_subspace_t tmp_st;
 double eigenvalue;
 diag_1d_t(propagator const& prop, sub_hilbert_space const& sp);
 void diag(double t, double dt);
};

// LAPACK diagonalization
struct diag_lapack_t {
 propagator const& prop;
 state_on_subspace_t from_st, to_st;
 matrix<dcomplex> workspace;
 std::pair<array<double, 1>, matrix<dcomplex>> eig;
 diag_lapack_t(propagator const& prop, sub_hilbert_space const& sp);
 void diag(double t, double dt);
};

class propagator {

 op_on_subspace_t const& h;        // Hamiltonian
 const bool is_static_h;           // Is Hamiltonian static on sp
 const h_interpolation h_interpol; // Hamiltonian interpolation
 const dcomplex h_coeff;           // Hamiltonian prefactor in the exponential
 const long N;                     // Dimension of Hilbert space

 // Choose the exact diagonaliztion solver
 enum {OneD, LAPACK, Lanczos} ed_solver;
 template<decltype(ed_solver) ed> using ed_type = std::integral_constant<decltype(ed_solver),ed>;

 friend class diag_1d_t;
 friend class diag_lapack_t;

 std::unique_ptr<diag_1d_t> diag_1d;
 std::unique_ptr<diag_lapack_t> diag_lapack;

 state_on_subspace_t apply_h(state_on_subspace_t const& st, double t, double dt) const;

public:

 using time_it_t = gf_mesh<retime>::const_iterator;

 propagator(op_on_subspace_t const& h, sub_hilbert_space const& sp,
            double hbar, bool is_static_op, h_interpolation h_interpol,
            long lanczos_min_matrix_size);

 void operator()(state_on_subspace_t & st,
                 time_it_t const& t_start, time_it_t const& t_end) const;

private:

 void propagate(state_on_subspace_t & st, time_it_t const& t_min, time_it_t const& t_max, dcomplex c) const;

 void propagate(state_on_subspace_t &, time_it_t const&, time_it_t const& t, dcomplex, ed_type<OneD>) const;
 void propagate(state_on_subspace_t &, time_it_t const&, time_it_t const& t, dcomplex, ed_type<LAPACK>) const;
 void propagate(state_on_subspace_t &, time_it_t const&, time_it_t const& t, dcomplex, ed_type<Lanczos>) const;
};

}
