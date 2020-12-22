/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <array>
#include <cstdlib>
#include <map>
#include <utility>

#include "types.hpp"
#include "init_state.hpp"
#include "worldlines.hpp"
#include "propagator.hpp"
#include "time_point_selector.hpp"

namespace realevol {

//
// Algorithm parameters for Lanczos propagation
//
struct lanczos_params_t {

  /// Use Lanczos algorithm to exponentiate matrices of this size or bigger
  int const& min_matrix_size;

  /// Lanczos convergence threshold for the energy, for each invariant subspace
  std::map<long, double> const& lanczos_energy_tol;

  /// Maximal dimension of the Krylov space, for each invariant subspace
  std::map<long, int> const& lanczos_max_krylov_dim;

  /// Lanczos energy tolerance on subspace 'spn'
  double get_energy_tol(int spn) const {
    auto it = lanczos_energy_tol.find(spn);
    return it != lanczos_energy_tol.end() ? it->second : -1;
  }

  /// Maximal dimension of the Krylov space within the 'spn' subspace
  int get_max_krylov_dim(int spn) const {
    auto it = lanczos_max_krylov_dim.find(spn);
    return it != lanczos_max_krylov_dim.end() ? it->second : -1;
  }
};

//
// Compute contributions to the dynamical trace of form
//
//  \Tr[\hat M_{N-1}(t_0) ... \hat M_1(t_{N-2}) \hat M_0(t_{N-1}) \hat\rho_0],
//
// where N = NPoints, M_n are Heisenberg operators taken at times t_n. Their
// corresponding Heisenberg operators must be monomials of creation/annihilation
// operators.
//
template<std::size_t NPoints, typename HamiltonianType>
class wl_worker {

  static_assert(NPoints != 0, "Trivial dynamical trace is not supported");

  using HScalarType = typename HamiltonianType::scalar_t;

  // Initial state of the system
  init_state const& initial_state;

  // Structure of the Hilbert space
  hilbert_space_structure<HamiltonianType> const& hss;

  // Are all matrix elements of H(t) time-independent on a given subspace?
  std::vector<bool> const& is_static_sp;

  // Selector of time point combinations to process.
  time_point_selector<NPoints> const& t_selector;

  // Hamiltonian H(t)
  op_on_subspace_t<HScalarType> h;

  // Planck's constant
  double hbar;

  // Time arguments of the operators (t_0, t_1, ..., t_{N-1}) in the trace are
  // defined on this mesh.
  gf_mesh<retime> t_mesh;

  // Hamiltonian interpolation method
  h_interpolation HInterpol;

  // Parameters for Lanczos propagation
  lanczos_params_t const& lanczos_params;

public:

  wl_worker(init_state const& initial_state,
            HamiltonianType const& h,
            double hbar,
            hilbert_space_structure<HamiltonianType> const& hss,
            std::vector<bool> const& is_static_sp,
            time_point_selector<NPoints> const& t_selector,
            gf_mesh<retime> const& t_mesh,
            h_interpolation HInterpol,
            lanczos_params_t const& lanczos_params);

  // Do the trace for one world line
  // The computed contribution will be *added* to 'result'
  void operator()(worldline_desc_t<NPoints> const& wl,
                  time_container_t<NPoints> & result) const;

private:

  using op_with_map_t = imperative_operator<sub_hilbert_space, dcomplex, true>;

  // Indices of time points of operators in the trace
  // Used to pass information between recursive calls to do_inner_loop<P>()
  mutable std::array<int, NPoints> t_indices;

  // Do an inner loop over a single time variable in the trace
  template<std::size_t Point>
  void do_inner_loop(std::array<propagator<HScalarType>, NPoints+1> const& props,
                     std::array<op_with_map_t, NPoints> const& ops,
                     std::array<state_on_subspace_t, NPoints+1> & psi,
                     state_on_subspace_t & phi,
                     dcomplex coeff,
                     time_container_t<NPoints> & result
                    ) const;
};

} // namespace realevol
