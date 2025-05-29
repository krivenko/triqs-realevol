// Copyright (c) 2013-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2013-2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2014-2016 Igor Krivenko
// Copyright (c) 2018 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Philipp D, Igor Krivenko, Olivier Parcollet, Nils Wentzell

#pragma once

#include <algorithm>
#include <complex>

#include <nda/nda.hpp>

#include <triqs/utility/exceptions.hpp>

namespace nda::lapack {

  template <bool Complex = false> struct tridiag_worker {};

  template <> struct tridiag_worker<false> {
    vector<double> D, E, W;
    matrix<double> Z;
    size_t s;
    tridiag_worker(size_t n) : s(0) { resize(n); }
    void resize(size_t n) {
      if (D.size() < n) {
        D.resize(n);
        E.resize(n);
        W.resize(std::max(1, 2 * int(n) - 2));
        Z.resize(n, n);
      }
      s = n;
    }
    template <typename VT> void operator()(VT const &diag, VT const &supdiag /* superdiagonal */) {
      if (supdiag.size() != diag.size() - 1) TRIQS_RUNTIME_ERROR << "Error in tridiagonal matrix diagonalization: size mismatch";
      resize(diag.size());
      std::copy(std::begin(diag), std::begin(diag) + s, std::begin(D));
      std::copy(std::begin(supdiag), std::begin(supdiag) + s - 1, std::begin(E));
      int info;
      f77::stev('V', s, D.data(), E.data(), Z.data(), first_dim(Z), W.data(), info);
      if (info != 0) TRIQS_RUNTIME_ERROR << "Error in tridiagonal matrix diagonalization " << info;
    }
    vector_view<double const> values() const {
      return vector_view<double const>({long(s)}, D.data());
    }
    matrix_view<double const> vectors() const {
      return matrix_view<double const>({long(s), long(s)}, Z.data());
    }
  };

  template <> struct tridiag_worker<true> {
    vector<double> D, E, W;
    matrix<double> Z;
    vector<std::complex<double>> U; // similarity transformation
    matrix<std::complex<double>> V; // transformed eigenvectors
    size_t s;
    tridiag_worker(size_t n) : s(0) { resize(n); }
    void resize(size_t n) {
      if (D.size() < n) {
        D.resize(n);
        E.resize(n);
        W.resize(std::max(1, 2 * int(n) - 2));
        Z.resize(n, n);
        U.resize(n);
        V.resize(n, n);
      }
      s = n;
    }
    template <typename VTd, typename VTe> void operator()(VTd const &diag, VTe const &supdiag /* superdiagonal */) {
      if (supdiag.size() != diag.size() - 1) TRIQS_RUNTIME_ERROR << "Error in tridiagonal matrix diagonalization: size mismatch";
      resize(diag.size());
      std::copy(std::begin(diag), std::begin(diag) + s, std::begin(D));
      // construct transformed off-diagonal elements;
      U(0) = 1.0;
      for (int i : range(0, s - 1)) {
        std::complex<double> e = supdiag(i);
        U(i + 1)               = U(i) * (std::isnormal(abs(e)) ? conj(e) / abs(e) : 1.0);
        E(i)                   = std::real(e * conj(U(i)) * (U(i + 1)));
      }
      int info;
      f77::stev('V', s, D.data(), E.data(), Z.data(), first_dim(Z), W.data(), info);
      if (info != 0) TRIQS_RUNTIME_ERROR << "Error in tridiagonal matrix diagonalization " << info;

      // Back-transform the eigenvectors
      for (int i : range(0, s))
        for (int j : range(0, s)) V(i, j) = Z(i, j) * conj(U(i)) * U(j);
    }
    vector_view<double const> values() const {
      return vector_view<double const>({long(s)}, D.data());
    }
    matrix_view<std::complex<double> const> vectors() const {
      return matrix_view<std::complex<double> const>({long(s), long(s)}, V.data());
    }
  };

} // namespace nda::lapack
