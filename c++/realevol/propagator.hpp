/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2022, I. Krivenko, M. Danilov, P. Kubiczek
 *
 * realevol is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * realevol is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * realevol. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

#include "hilbert_space/hilbert_space.hpp"
#include "hilbert_space/state.hpp"
#include "hilbert_space/imperative_operator.hpp"

#include "types.hpp"
#include "lanczos_worker.hpp"

namespace realevol {

// Hamiltonian interpolation between time slices
enum h_interpolation {Rectangle, Trapezoid, Simpson};

template<typename HScalarType> class propagator;

// Instantaneous Hamiltonian
template<typename HScalarType>
struct inst_h_t {
 const bool is_static;
 op_on_subspace_t<HScalarType> const& h;        // Hamiltonian
 const h_interpolation h_interpol; // Hamiltonian interpolation
 double t;
 double dt;
 inst_h_t(op_on_subspace_t<HScalarType> const& h, bool is_static, h_interpolation h_interpol);
 inst_h_t(inst_h_t const&) = delete;
 inst_h_t(inst_h_t &&) noexcept = default;
 state_on_subspace_t operator()(state_on_subspace_t const& st) const;
};

// Propagation on 1-d subspaces
template<typename HScalarType>
struct prop_1d_t {
 inst_h_t<HScalarType> inst_h;
 state_on_subspace_t tmp_st;
 double eigenvalue = {};
 prop_1d_t(inst_h_t<HScalarType> && inst_h, sub_hilbert_space const& sp);
 prop_1d_t(prop_1d_t const&) = delete;
 void diag(double t, double dt);
 void operator()(state_on_subspace_t & st, double t_min, double t_max, dcomplex c);
};

// LAPACK propagation
template<typename HScalarType>
struct prop_lapack_t {
 inst_h_t<HScalarType> inst_h;
 state_on_subspace_t from_st, to_st;
 matrix<dcomplex> workspace;
 std::pair<array<double, 1>, matrix<dcomplex>> eig;
 prop_lapack_t(inst_h_t<HScalarType> && inst_h, sub_hilbert_space const& sp);
 prop_lapack_t(prop_lapack_t const&) = delete;
 void diag(double t, double dt);
 void operator()(state_on_subspace_t & st, double t_min, double t_max, dcomplex c);
};

// Lanczos propagation
template<typename HScalarType>
struct prop_lanczos_t {
 inst_h_t<HScalarType> inst_h;
 lanczos_worker<inst_h_t<HScalarType>, state_on_subspace_t> lw;
 vector<dcomplex> krylov_coeffs;
 prop_lanczos_t(inst_h_t<HScalarType> && inst_h, sub_hilbert_space const& sp,
                double gs_energy_convergence, int max_krylov_dim);
 prop_lanczos_t(prop_lanczos_t const&) = delete;
 void operator()(state_on_subspace_t & st, double t_min, double t_max, dcomplex c);
};

template<typename HScalarType>
class propagator {

 const bool is_static_h;           // Is Hamiltonian static on sp
 const dcomplex h_coeff;           // Hamiltonian prefactor in the exponential
 gf_mesh<retime> const& t_mesh;    // Time mesh

 // Choose the exact diagonaliztion solver
 enum {OneD, LAPACK, Lanczos} ed_solver;

 std::shared_ptr<void> prop_impl;

 void propagate(state_on_subspace_t & st, int t_min_index, int t_max_index, dcomplex c) const;

public:

 propagator(op_on_subspace_t<HScalarType> const& h, sub_hilbert_space const& sp,
            gf_mesh<retime> const& t_mesh,
            double hbar, bool is_static_op, h_interpolation h_interpol,
            long lanczos_min_matrix_size, double gs_energy_convergence, int max_krylov_dim);

 void operator()(state_on_subspace_t & st, int t_start_index, int t_end_index) const;
};

class time_expr;
class time_interp;

extern template class propagator<time_expr>;
extern template class propagator<time_interp>;

}
