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

namespace realevol {

// Replace a callable object with a constant if it returns equal values at all mesh points
template<class Callable, class Mesh>
Callable try_reduce_to_constant(Callable const& callable, Mesh const& m) {
  auto it = m.begin();
  auto value = callable(*it);

  using triqs::utility::is_zero;
  for(it++; it != m.end(); it++)
    if(!is_zero(callable(*it) - value)) return callable;

  return is_zero(value.imag()) ? Callable(value.real()) : Callable(value);
}

// Predicate to check whether a callable object returns zero on a mesh
template<class Mesh>
struct is_zero_on_mesh {

 Mesh const& mesh;
 is_zero_on_mesh(Mesh const& mesh) : mesh(mesh) {}

 template<class Callable> bool operator()(Callable const& te) {
  using triqs::utility::is_zero;
  if(te.is_zero()) return true;
  for(auto t : mesh)
   if(!is_zero(te(t))) return false;
  return true;
 }
};

}
