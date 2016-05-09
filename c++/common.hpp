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

#include <complex>
#include <type_traits>

#include "time_expr.hpp"

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/space_partition.hpp>
#include <triqs/utility/numeric_ops.hpp>

namespace realevol {

using realevol::operators::many_body_operator_generic;
using realevol::hilbert_space::imperative_operator;
using realevol::hilbert_space::fock_state_t;
using realevol::hilbert_space::state;
using realevol::hilbert_space::fundamental_operator_set;
using realevol::hilbert_space::hilbert_space;
using realevol::hilbert_space::sub_hilbert_space;
using realevol::hilbert_space::space_partition;
using realevol::hilbert_space::project;

using operator_t = many_body_operator_generic<time_expr>;
using hilbert_space_t = hilbert_space::hilbert_space;

using op_on_space_t = imperative_operator<hilbert_space_t,time_expr,false>;
using op_on_subspace_t = imperative_operator<sub_hilbert_space,time_expr,false>;
using state_on_space_t = state<hilbert_space_t,dcomplex,true>;
using state_on_subspace_t = state<sub_hilbert_space,dcomplex,false>;
using dyn_state_on_space_t = state<hilbert_space_t,time_expr,true>;

}
