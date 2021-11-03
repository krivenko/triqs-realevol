/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2021, I. Krivenko, M. Danilov, P. Kubiczek
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
#include <utility>

namespace realevol {

/// Apply a given functor 'l' to an integer sequence 0, 1, ... N-1 and return
/// a new std::array made of the return values of 'l'
template<typename T, typename Lambda, std::size_t... Ints>
std::array<T, sizeof...(Ints)> make_array_impl(Lambda &&l, std::index_sequence<Ints...>) {
  return {l(Ints)...};
}
template<typename T, std::size_t N, typename Lambda>
std::array<T, N> make_array(Lambda &&l) {
  return make_array_impl<T>(std::forward<Lambda>(l), std::make_index_sequence<N>());
}

/// Apply a given functor 'l' to all elements of std::array 'a' and return
/// a new std::array made of the return values of 'l'
template<typename Out, typename In, typename Lambda, std::size_t... Ints>
std::array<Out, sizeof...(Ints)> map_array_impl(Lambda &&l,
                                                std::array<In, sizeof...(Ints)> const &a,
                                                std::index_sequence<Ints...>) {
  return {l(a[Ints])...};
}
template<typename Out, typename In, std::size_t N, typename Lambda>
std::array<Out, N> map_array(Lambda &&l, std::array<In, N> const &a) {
  return map_array_impl<Out>(std::forward<Lambda>(l), a, std::make_index_sequence<N>());
}

/// Element-wise composition of arrays 'a' and 'b' using a given functor 'l'
template<typename Out, typename In1, typename In2, typename Lambda, std::size_t... Ints>
std::array<Out, sizeof...(Ints)>
map_array_impl2(Lambda &&l, std::array<In1, sizeof...(Ints)> const &a,
                std::array<In2, sizeof...(Ints)> const &b,
                std::index_sequence<Ints...>) {
  return {l(a[Ints], b[Ints])...};
}
template<typename Out, typename In1, typename In2, std::size_t N, typename Lambda>
std::array<Out, N> map_array2(Lambda &&l, std::array<In1, N> const &a, std::array<In2, N> const &b) {
  return map_array_impl2<Out>(std::forward<Lambda>(l), a, b, std::make_index_sequence<N>());
}

/// Make a new std::array<T, N> with all elements equal to T(args...)
template<typename T, std::size_t N, typename... Args> std::array<T, N>
make_array_repeat(Args &&... args) {
  return make_array<T, N>([&](std::size_t) { return T(args...); });
}

} // namespace realevol
