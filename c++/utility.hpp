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

#include <utility>

#include "common.hpp"

namespace realevol {

/// Compute retarded and advanced Green's functions out of the lesser and greater components
/**
   @param g_l Lesser Green's function
   @param g_l Greater Green's function
   @return Retarded and advanced Green's functions
 */
std::pair<block_gf_2t_t,block_gf_2t_t> make_gf_ret_adv(block_gf_2t_t const& g_l, block_gf_2t_t const& g_g);

}

