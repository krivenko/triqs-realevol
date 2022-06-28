/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2023, I. Krivenko, M. Danilov, P. Kubiczek
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
#include <cmath>
#include <utility>

#include <triqs/gfs.hpp>
#include <triqs/mesh/retime.hpp>

namespace realevol {

/// Select (t_0, t_1, t_{N-1}) tuples fulfilling some inequalities
template<std::size_t NPoints> class time_point_selector {

protected:

 const std::array<std::pair<double, double>, NPoints> t_ranges;
 const std::array<double, NPoints - 1> delta_t_max;

 using mesh_point_t = triqs::mesh::retime::mesh_point_t;

 [[nodiscard]] inline bool in_domain(std::array<mesh_point_t, NPoints> const& t) const {
   // Check the ranges first
   for(std::size_t p = 0; p < NPoints; ++p) {
     if(t[p] < t_ranges[p].first || t[p] > t_ranges[p].second)
       return false;
   }
   // Then check the time separations
   for(std::size_t p = 0; p < NPoints - 1; ++p) {
     if(std::abs(double(t[p]) - double(t[p + 1])) > delta_t_max[p])
       return false;
   }
   return true;
 }

public:

 time_point_selector(std::array<std::pair<double, double>, NPoints> t_ranges,
                     std::array<double, NPoints - 1> delta_t_max
                    ) :
  t_ranges(std::move(t_ranges)),
  delta_t_max(std::move(delta_t_max)) {}

 inline bool operator()(std::array<mesh_point_t, NPoints> const& t) const {
    return in_domain(t);
  }
};

/// Select (t_0, t_1) tuples fulfilling some inequalities and the t_0 >= t_1 condition
class time_point_selector_lower_triangle : public time_point_selector<2> {

public:

  time_point_selector_lower_triangle(
    std::array<std::pair<double, double>, 2> t_ranges,
    std::array<double, 1> delta_t_max
  ) :
  time_point_selector<2>(std::move(t_ranges), std::move(delta_t_max)) {}

  inline bool operator()(std::array<mesh_point_t, 2> const& t) const {
    return time_point_selector<2>::operator()(t) && (t[0] >= t[1]);
  }
};

}
