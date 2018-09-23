/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018 I. Krivenko
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

#include <algorithm>
#include <cmath>
#include <vector>

#include <triqs/gfs/meshes/segment.hpp>
#include <triqs/utility/exceptions.hpp>

namespace realevol {

template<class T>
class interpolator1d {
  public:
  using value_type = T;
  using mesh_type = triqs::gfs::segment_mesh;
  using data_type = std::vector<value_type>;

private:

  mesh_type const& mesh;
  data_type data;
  bool is_constant_;

public:

  interpolator1d() = default;
  interpolator1d(mesh_type const& mesh, data_type const& data) : mesh(mesh), data(data) {
    if(mesh.size() != data.size()) TRIQS_RUNTIME_ERROR << "Inconsistent sizes of mesh and data";
    is_constant_ = std::adjacent_find(data.begin(), data.end(), std::not_equal_to<value_type>()) == data.end();
  }
  interpolator1d(mesh_type const& mesh, value_type val) : mesh(mesh), data(mesh.size(), val), is_constant_(true) {}

  value_type operator()(double x) const {

    TRIQS_ASSERT(x >= mesh.x_min());
    TRIQS_ASSERT(x <= mesh.x_max());

    if(is_constant_) return data[0];

    int l = std::floor(x / mesh.delta());
    double d = (x - mesh[l]) / (mesh[l + 1] - mesh[l]);
    return (1.0 - d) * data[l] + d * data[l + 1];
  }

  // Accessors
  mesh_type const& get_mesh() const { return mesh; }
  data_type const& get_data() const { return data; }

  friend bool is_constant(interpolator1d const& interp) { return interp.is_constant_; }
  bool is_zero() const { return is_constant(*this) && data[0] == T{}; }
};

}
