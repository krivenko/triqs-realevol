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
#include <triqs/gfs.hpp>

namespace realevol {

/// Select (t,tp) pairs fulfilling some inequalities
class time_point_selector {

 const std::pair<double, double> t_range;
 const std::pair<double, double> tp_range;
 const double delta_t_max;

 using mesh_point_t = triqs::gfs::gf_mesh<triqs::gfs::retime>::mesh_point_t;

 inline bool in_domain(mesh_point_t const &t, mesh_point_t const &tp) const {
   return (t >= t_range.first && t <= t_range.second) &&
          (tp >= tp_range.first && tp <= tp_range.second) &&
          std::abs(double(t) - double(tp)) <= delta_t_max;
 }

public:

 time_point_selector(std::pair<double, double> const& t_range,
                     std::pair<double, double> const& tp_range,
                     double delta_t_max) :
  t_range(t_range), tp_range(tp_range), delta_t_max(delta_t_max) {}

 inline bool operator()(mesh_point_t const &t, mesh_point_t const &tp) const {
  bool b = in_domain(t, tp);
#ifdef USE_ANTIHERMITICITY
  if(t >= tp) { // lower triangle: check that the (t,tp) pair is within domain
   return b;
  } else { // upper triangle: check that (t,tp) does not have an equivalent pair in the lower triangle
   return b && !in_domain(tp, t);
  }
#else
  return b;
#endif
 }
};

}
