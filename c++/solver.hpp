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

#include <string>
#include <vector>
#include <map>
#include <utility>

#include <triqs/mpi/base.hpp>
#include <triqs/gfs.hpp>

#include "time_expr.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include "solve_parameters.hpp"

namespace realevol {

using namespace triqs::gfs;

using indices_type = triqs::operators::indices_t;
using g_2t_t = gf<cartesian_product<retime, retime>>;

class solver {

 std::map<std::string, indices_type> gf_struct;           // Block structure of the Green function
 block_gf<cartesian_product<retime,retime>> g_ret, g_adv; // Advanced and retarded GF's to be calculated
 triqs::mpi::communicator comm;                           // MPI communicator
 solve_parameters_t params;                               // Parameters of the last call to solve

public:

 solver(std::map<std::string,indices_type> const& gf_struct, std::pair<double,double> time_window, int n_t = 1000);

 TRIQS_WRAP_ARG_AS_DICT
 /// Solve the many-bosy problem for the given Hamiltonian h
 void solve(solve_parameters_t const& p);

 /// Set of parameters used in the last call to solve
 solve_parameters_t get_last_run_parameters() const { return params; }

 /// Retarded GF in real time
 block_gf_view<cartesian_product<retime,retime>> get_g_ret() { return g_ret; }

 /// Advanced GF in real time
 block_gf_view<cartesian_product<retime,retime>> get_g_adv() { return g_adv; }
};

}
