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

#include <utility>

#include "types.hpp"

namespace realevol {

/// Compute retarded and advanced Green's functions out of the lesser and greater components
/**
   @param g_l Lesser Green's function
   @param g_g Greater Green's function
   @return Retarded and advanced Green's functions
 */
std::pair<block_gf_2t_t,block_gf_2t_t> make_gf_ret_adv(block_gf_2t_t const& g_l,
                                                       block_gf_2t_t const& g_g);

}

