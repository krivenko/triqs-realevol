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

#include <array>
#include <cmath>
#include <utility>

#include <triqs/gfs.hpp>

namespace realevol {

/// Select (t_0, t_1, t_{N-1}) tuples fulfilling some inequalities
template<std::size_t NPoints> class time_point_selector {

 const std::array<std::pair<double, double>, NPoints> t_ranges;
 const std::array<double, NPoints - 1> delta_t_max;
 // TODO: Remove these if possible
 bool use_antihermiticity;
 bool swap_t_tp;

 using mesh_point_t = triqs::gfs::gf_mesh<triqs::gfs::retime>::mesh_point_t;

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
                     std::array<double, NPoints - 1> delta_t_max,
                     bool use_antihermiticity = false,
                     bool swap_t_tp = false
                    ) :
  t_ranges(std::move(t_ranges)),
  delta_t_max(std::move(delta_t_max)),
  use_antihermiticity(use_antihermiticity),
  swap_t_tp(swap_t_tp) {}

 void set_swap_t_tp(bool swap) { swap_t_tp = swap; }

 inline bool operator()(std::array<mesh_point_t, NPoints> const& t) const {
    if constexpr(NPoints == 2) {
      auto t1 = swap_t_tp ? t[1] : t[0];
      auto t2 = swap_t_tp ? t[0] : t[1];
      bool b = in_domain({t1, t2});
      if(use_antihermiticity) {
        if(t1 >= t2) { // lower triangle: check that the (t,tp) pair is within domain
          return b;
        } else { // upper triangle: check that (t,tp) does not have an equivalent pair in the lower triangle
          return b && !in_domain({t2, t1});
        }
      } else
        return in_domain(t);
    } else {
      return in_domain(t);
    }
  }
};

}
