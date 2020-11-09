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
#include <iostream>
#include <vector>

#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>
#include <triqs/utility/exceptions.hpp>

namespace realevol {

template<class T>
class interpolator1d {

public:
  using value_type = T;
  using mesh_type = triqs::gfs::segment_mesh;
  using data_type = triqs::arrays::array<value_type, 1>;

protected:

  bool is_constant_;
  mesh_type mesh;
  data_type data;

public:

  // Construct explicitly constant interpolator
  interpolator1d() : is_constant_(true), data{0} {};
  interpolator1d(value_type val) : is_constant_(true), data{val} {}
  interpolator1d(mesh_type const& mesh, value_type val) :
    is_constant_(true),
    mesh(mesh),
    data{val}
  {}

  // Construct an interpolator from a mesh and an array of values
  interpolator1d(mesh_type const& mesh, data_type const& data_) :
    is_constant_(std::adjacent_find(data_.begin(), data_.end(),
                                    std::not_equal_to<value_type>()) == data_.end()),
    mesh(mesh),
    data(is_constant_ ? data_type{data_(0)} : data_) {
    // Optionally save memory by storing only one data value
    if(!is_constant_ && mesh.size() != data_.shape()[0])
      TRIQS_RUNTIME_ERROR << "Inconsistent sizes of mesh and data";
    if(mesh.size() < 2)
      TRIQS_RUNTIME_ERROR << "Too few data points, need at least 2";
  }

  value_type operator()(double x) const {

    if(is_constant_) return data(0);

    TRIQS_ASSERT(x >= mesh.x_min());
    TRIQS_ASSERT(x <= mesh.x_max());

    int l = std::floor(x / mesh.delta());
    if(l == mesh.size() - 1) return data(l);
    double d = (x - mesh[l]) / (mesh[l + 1] - mesh[l]);
    return (1.0 - d) * data(l) + d * data(l + 1);
  }

  size_t size() const { return mesh.size(); }

  // Accessors
  mesh_type const& get_mesh() const { return mesh; }
  data_type const& get_data() const { return data; }
  data_type & get_data() { return data; }

  bool operator==(interpolator1d const& interp) const {
    if(is_constant_) {
      return interp.is_constant_ && data(0) == interp.data(0);
    } else {
      return (!interp.is_constant_) && mesh == interp.mesh && data == interp.data;
    }
  }

  friend bool is_constant(interpolator1d const& interp) { return interp.is_constant_; }
  bool is_zero() const { return is_constant(*this) && data(0) == T{}; }

  friend std::ostream& operator<<(std::ostream& os, interpolator1d const& interp) {
    if(interp.is_constant_)
      os << interp.data(0);
    else {
      os << "ti([" << interp.mesh.x_min() << "," << interp.mesh.x_max() << "]->[";
      os << interp.data(0) << (interp.mesh.size() > 2 ? ",...," : ",")
         << interp.data(interp.mesh.size()-1) << "])";
    }
    return os;
  }

};

}
