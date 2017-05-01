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
#include <functional>

#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include <triqs/utility/variant.hpp>

#include "common.hpp"
#include "lanczos_worker.hpp"

namespace realevol {

// Hamiltonian interpolation between time slices
enum h_interpolation {Rectangle, Trapezoid, Simpson};

using time_it_t = gf_mesh<retime>::const_iterator;

class propagator;

// Instantaneous Hamiltonian
struct inst_h_t {
 const bool is_static;
 op_on_subspace_t const& h;        // Hamiltonian
 const h_interpolation h_interpol; // Hamiltonian interpolation
 double t;
 double dt;
 inst_h_t(op_on_subspace_t const& h, bool is_static, h_interpolation h_interpol);
 inst_h_t(inst_h_t const&) = delete;
 inst_h_t(inst_h_t &&) noexcept = default;
 state_on_subspace_t operator()(state_on_subspace_t const& st) const;
};

// Propagation on 1-d subspaces
struct prop_1d_t {
 inst_h_t inst_h;
 state_on_subspace_t tmp_st;
 double eigenvalue;
 prop_1d_t(inst_h_t && inst_h, sub_hilbert_space const& sp);
 prop_1d_t(prop_1d_t const&) = delete;
 void diag(double t, double dt);
 void operator()(state_on_subspace_t & st,
                 time_it_t const& t_min, time_it_t const& t_max,
                 dcomplex c);
};

// LAPACK propagation
struct prop_lapack_t {
 inst_h_t inst_h;
 state_on_subspace_t from_st, to_st;
 matrix<dcomplex> workspace;
 std::pair<array<double, 1>, matrix<dcomplex>> eig;
 prop_lapack_t(inst_h_t && inst_h, sub_hilbert_space const& sp);
 prop_lapack_t(prop_lapack_t const&) = delete;
 void diag(double t, double dt);
 void operator()(state_on_subspace_t & st,
                 time_it_t const& t_min, time_it_t const& t_max,
                 dcomplex c);
};

// Lanczos propagation
struct prop_lanczos_t {
 inst_h_t inst_h;
 lanczos_worker<inst_h_t, state_on_subspace_t> lw;
 matrix<dcomplex> lanczos_exp;
 prop_lanczos_t(inst_h_t && inst_h, sub_hilbert_space const& sp,
                double gs_energy_convergence, int max_krylov_dim);
 prop_lanczos_t(prop_lanczos_t const&) = delete;
 void operator()(state_on_subspace_t & st,
                 time_it_t const& t_min, time_it_t const& t_max,
                 dcomplex c);
};

class propagator {

 const bool is_static_h;           // Is Hamiltonian static on sp
 const dcomplex h_coeff;           // Hamiltonian prefactor in the exponential

 // Choose the exact diagonaliztion solver
 enum {OneD, LAPACK, Lanczos} ed_solver;

 std::shared_ptr<void> prop_impl;
 void propagate(state_on_subspace_t & st, time_it_t const& t_min, time_it_t const& t_max, dcomplex c) const;

public:

 propagator(op_on_subspace_t const& h, sub_hilbert_space const& sp,
            double hbar, bool is_static_op, h_interpolation h_interpol,
            long lanczos_min_matrix_size, double gs_energy_convergence, int max_krylov_dim);

 void operator()(state_on_subspace_t & st,
                 time_it_t const& t_start, time_it_t const& t_end) const;
};

}
