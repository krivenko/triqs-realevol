/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2026, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include <triqs/gfs.hpp>
#include <triqs/mesh/retime.hpp>
#include <triqs/mesh/prod.hpp>
#include <triqs/utility/numeric_ops.hpp>

#include "operators/many_body_operator.hpp"
#include "hilbert_space/fundamental_operator_set.hpp"
#include "hilbert_space/hilbert_space.hpp"
#include "hilbert_space/imperative_operator.hpp"
#include "hilbert_space/state.hpp"

namespace realevol {

using namespace triqs;
using namespace triqs::gfs;
using namespace realevol::operators;
using namespace realevol::hilbert_space;

namespace detail {

template<typename T, std::size_t...> using return_T = T;

template<typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
struct make_cartesian_product;

template<typename T, std::size_t N, std::size_t... Indices>
struct make_cartesian_product<T, N, std::index_sequence<Indices...>> {
  using type = mesh::prod<return_T<T, Indices>...>;
};

template<typename T> struct make_cartesian_product<T, 1> {
  using type = T;
};

} // namespace detail

// Time meshes
template<std::size_t NPoints>
using time_mesh_t = typename detail::make_cartesian_product<mesh::retime, NPoints>::type;
using mesh_t_t = time_mesh_t<1>;
using mesh_2t_t = time_mesh_t<2>;
using mesh_3t_t = time_mesh_t<3>;

// Multi-dimensional GF containers for computed correlation functions
template<std::size_t NPoints>
using time_container_t = gf<typename detail::make_cartesian_product<mesh::retime, NPoints>::type,
                            scalar_valued>;
using expectval_container_t = time_container_t<1>;
using correlator_2t_container_t = time_container_t<2>;
using correlator_3t_container_t = time_container_t<3>;

using gf_2t_t = gf<mesh::prod<mesh::retime, mesh::retime>>;
using block_gf_2t_t = block_gf<mesh::prod<mesh::retime, mesh::retime>>;
using gf_2t_view = gf_view<mesh::prod<mesh::retime, mesh::retime>>;
using block_gf_2t_view = block_gf_view<mesh::prod<mesh::retime, mesh::retime>>;

using indices_type = operators::indices_t;
using chi_indices_t = std::vector<std::pair<std::string, long>>;

template<std::size_t NPoints>
using time_container_view_t = typename time_container_t<NPoints>::view_type;

using static_operator_t = many_body_operator;

template<typename ScalarType>
using op_on_space_t = imperative_operator<class hilbert_space, ScalarType, false>;
template<typename ScalarType>
using op_on_subspace_t = imperative_operator<sub_hilbert_space, ScalarType, false>;

using static_op_on_space_t = imperative_operator<class hilbert_space, dcomplex, false>;
using static_op_on_subspace_t = imperative_operator<sub_hilbert_space, dcomplex, false>;

using state_on_space_t = state<class hilbert_space, dcomplex, true>;
using state_on_subspace_t = state<sub_hilbert_space, dcomplex, false>;
template<typename ScalarType>
using dyn_state_on_space_t = state<class hilbert_space, ScalarType, true>;

} // namespace realevol
