/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <triqs/gfs.hpp>
#include <triqs/utility/numeric_ops.hpp>

#include "operators/many_body_operator.hpp"
#include "hilbert_space/fundamental_operator_set.hpp"
#include "hilbert_space/hilbert_space.hpp"
#include "hilbert_space/imperative_operator.hpp"
#include "hilbert_space/state.hpp"

#include "time_expr.hpp"
#include "time_interp.hpp"

namespace realevol {

using namespace triqs::gfs;
using namespace realevol::operators;
using namespace realevol::hilbert_space;

using gf_2t_t = gf<cartesian_product<retime, retime>>;
using block_gf_2t_t = block_gf<cartesian_product<retime, retime>>;
using gf_2t_view = gf_view<cartesian_product<retime, retime>>;
using block_gf_2t_view = block_gf_view<cartesian_product<retime, retime>>;

using time_expr_operator_t = realevol::operators::many_body_operator_generic<time_expr>;
using time_interp_operator_t = realevol::operators::many_body_operator_generic<time_interp>;
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

}
