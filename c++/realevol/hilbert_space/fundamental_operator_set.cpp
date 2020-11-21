/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013-2020, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include "hdf5_hl.h"

// FIXME: Code in <triqs/utility/variant_extensions.hpp> depends on these
// headers but does not include them.
//
// https://github.com/TRIQS/triqs/pull/799
#include <ostream>
#include <sstream>
#include <vector>

#include <triqs/utility/variant_extensions.hpp>
#include <h5/h5.hpp>

#include "fundamental_operator_set.hpp"

namespace realevol::hilbert_space {


 namespace { // auxiliary functions

  // a little visitor for reduction to string
  struct variant_visitor {
   std::string operator()(int i) const { return "i" + std::to_string(i); }
   std::string operator()(std::string const& s) const { return "s" + s; }
  };

  // decode the string
  std::variant<int, std::string> string_to_variant(std::string const& s) {
   switch (s[0]) {
    case 'i':
     return std::stoi(s.c_str() + 1); // the variant is an int. Skip the first char and recover the int
    case 's':
     return s.c_str() + 1; // the variant is a string. Just skip the first char
    default:
     TRIQS_RUNTIME_ERROR << "Variant indices absent in h5 read";
   }
  }

  // fundamental_operator_set --> vec vec string
  std::vector<std::vector<std::string>> to_vec_vec_string(fundamental_operator_set const& fops, statistic_enum stat) {
   std::vector<std::vector<std::string>> v(fops.size(stat));
   for (auto it = fops.begin(stat); it != fops.end(stat); ++it) { // loop over the couple (indices list, number)
    if (it->linear_index >= fops.size()) TRIQS_RUNTIME_ERROR << " Internal error fundamental_operator_set to vec vec string";
    for (auto& x : it->index) v[it->linear_index].push_back(visit(variant_visitor{}, x));
    // variants x are transformed to a string, add 'i' or 's' in front of the string
   }
   return v;
  }

  fundamental_operator_set::indices_t to_indices(std::vector<std::string> const& v) {
   fundamental_operator_set::indices_t indices; // list of indices of this C, C^+ op
   for (auto& x : v)
    if (!x.empty()) indices.push_back(string_to_variant(x));
   return indices;
  }

 } // auxiliary functions


 // private constructor
 fundamental_operator_set::fundamental_operator_set(std::vector<std::vector<std::string>> const& v_f,
                                                    std::vector<std::vector<std::string>> const& v_b) {
  for (int n = 0; n < v_f.size(); ++n) fermion_map_index_n.insert({to_indices(v_f[n]), n});
  for (int n = 0; n < v_b.size(); ++n) boson_map_index_n.insert({to_indices(v_b[n]), n});
 }

 // --- h5

 void h5_write_attribute(h5::object obj, std::string const& name, fundamental_operator_set const& fops) {
  h5::h5_write_attribute(obj, name, to_vec_vec_string(fops, Fermion));
  if(fops.size(Boson) != 0)
   h5::h5_write_attribute(obj, name + "_boson", to_vec_vec_string(fops, Boson));
 }

 void h5_read_attribute(h5::object obj, std::string const& name, fundamental_operator_set & fops) {
  std::vector<std::vector<std::string>> v_f;
  h5::h5_read_attribute(obj, name, v_f);
  std::string name_boson = name + "_boson";
  if(H5LTfind_attribute(obj, name_boson.c_str()) != 0) {
   std::vector<std::vector<std::string>> v_b;
   h5::h5_read_attribute(obj, name_boson, v_b);
   fops = fundamental_operator_set(v_f, v_b);
  } else
   fops = fundamental_operator_set(v_f);
 }

}
