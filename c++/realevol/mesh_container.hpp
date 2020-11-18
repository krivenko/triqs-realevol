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

#include <vector>
#include <ostream>
#include <type_traits>

#include <boost/iterator/iterator_facade.hpp>

#include <triqs/utility/exceptions.hpp>

namespace realevol {

template<class T, class Mesh>
class mesh_container : public std::vector<T> {

public:

 using mesh_t = Mesh;
 using base_t = std::vector<T>;
 using value_type = typename base_t::value_type;

 // vector-compatible constructors
 template<class... ValueConstructorArgs>
 explicit mesh_container(mesh_t const& mesh, ValueConstructorArgs && ...vc_args) :
  base_t(mesh.size(),value_type(vc_args...)), mesh(mesh) {}

 mesh_container(mesh_t const& mesh, base_t const& value) : base_t(value), mesh(mesh) {}

 mesh_t const& get_mesh() const { return mesh; }

private:

 template<bool Const> struct point_value_pair {
  typename mesh_t::mesh_point_t mesh_point;
  std::conditional_t<Const,value_type const&, value_type&> value;
 };

 template<bool Const> struct iterator_ : public boost::iterator_facade<iterator_<Const>,
 point_value_pair<Const>, boost::random_access_traversal_tag, point_value_pair<Const>> {
  typename mesh_t::const_iterator m_it;
  std::conditional_t<Const,typename base_t::const_iterator,typename base_t::iterator> v_it;
 public :
  iterator_(decltype(m_it) m_it, decltype(v_it) v_it) : m_it(m_it), v_it(v_it) {}
  iterator_() {}
  iterator_(iterator_ const& it) = default;
  iterator_(iterator_ && it) = default;
  iterator_ & operator=(iterator_ const& it) = default;
  iterator_ & operator=(iterator_ && it) noexcept {
   using std::swap;
   swap(it.m_it, this->m_it);
   swap(it.v_it, this->v_it);
   return *this;
  }
 private:
  friend class boost::iterator_core_access;
  point_value_pair<Const> dereference() const { return {*m_it,*v_it}; }
  bool equal(iterator_ const& it) const { return m_it == it.m_it && v_it == it.v_it; }
  void increment() { ++m_it; ++v_it; }
  void decrement() { --m_it; --v_it; }
  void advance(std::size_t n) { m_it += n; v_it += n; }
  std::size_t distance_to(iterator_ const& it) const { return std::distance(it.m_it, m_it); }
 };

public:

 // Iterator over pairs mesh point-value pairs (const)
 using const_iterator = iterator_<true>;
 const_iterator begin() const noexcept {
  return const_iterator(std::begin(mesh), std::begin(static_cast<base_t const&>(*this))); }
 const_iterator cbegin() const noexcept {
  return const_iterator(std::begin(mesh), std::begin(static_cast<base_t const&>(*this))); }
 const_iterator end() const noexcept {
  return const_iterator(std::end(mesh), std::end(static_cast<base_t const&>(*this))); }
 const_iterator cend() const noexcept {
  return const_iterator(std::end(mesh), std::end(static_cast<base_t const&>(*this))); }

 // Iterator over pairs mesh point-value pairs
 using iterator = iterator_<false>;
 iterator begin() noexcept { return iterator(std::begin(mesh), std::begin(static_cast<base_t&>(*this))); }
 iterator end() noexcept { return iterator(std::end(mesh), std::end(static_cast<base_t&>(*this))); }

 // Insert contents of the container into a stream as two columns
 friend std::ostream& operator<<(std::ostream & os, mesh_container & MC) {
  for(auto e : MC) os << e.mesh_point << '\t' << e.value << std::endl;
  return os;
 }

private:
 mesh_t mesh;

 /// Write into HDF5
 friend void h5_write(h5::group fg, std::string subgroup_name, mesh_container const &c) {
  h5::group gr = fg.create_group(subgroup_name);
  h5_write(gr, "vector", static_cast<base_t const&>(c));
  h5_write(gr, "mesh", c.mesh);
 }
 /// Read from HDF5
 friend void h5_read(h5::group fg, std::string subgroup_name, mesh_container &c) {
  h5::group gr = fg.open_group(subgroup_name);
  h5_read(gr, "vector", static_cast<base_t&>(c));
  h5_read(gr, "mesh", c.mesh);
  if(c.size() != c.mesh.size())
   TRIQS_RUNTIME_ERROR << "Error reading from HDF5! Inconsistent sizes of the mesh and the data.";
 }

};

}
