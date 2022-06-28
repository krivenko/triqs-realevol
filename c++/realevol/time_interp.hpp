/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2022, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <limits>

#include <boost/operators.hpp>

#include "interpolator.hpp"

namespace realevol {

class time_interp : public interpolator1d<std::complex<double>>,
                    public boost::field_operators<time_interp> {

public:

  using interp_t = interpolator1d<std::complex<double>>;
  using real_data_type = triqs::arrays::array<double, 1>;

  // Constructors
  time_interp() = default;
  time_interp(double r);
  time_interp(double r, double i);
  time_interp(std::complex<double> const& z);

  time_interp(mesh_type const& mesh);
  time_interp(mesh_type const& mesh, double r);
  time_interp(mesh_type const& mesh, double r, double i);
  time_interp(mesh_type const& mesh, std::complex<double> const& z);
  time_interp(mesh_type const& mesh, real_data_type const& r);
  time_interp(mesh_type const& mesh, real_data_type const& r, double i);
  time_interp(mesh_type const& mesh, double r, real_data_type const& i);
  time_interp(mesh_type const& mesh, real_data_type const& r, real_data_type const& i);
  time_interp(mesh_type const& mesh, data_type const& z);

  time_interp(time_interp const&) = default;
  time_interp(time_interp &&) noexcept = default;
  ~time_interp() = default;

  // Assignment
  time_interp & operator=(time_interp const& ti) = default;
  time_interp & operator=(time_interp && ti) noexcept = default;

  // Arithmetics
  time_interp operator-() const;

  time_interp & operator+=(time_interp const& ti);
  time_interp & operator-=(time_interp const& ti);
  time_interp & operator*=(time_interp const& ti);
  time_interp & operator/=(time_interp const& ti);

  // Mixed arithmetics time_interp + double
  time_interp & operator+=(double r);
  friend time_interp operator+(time_interp const& ti, double r);
  friend time_interp operator+(double r, time_interp const& ti);

  time_interp & operator-=(double r);
  friend time_interp operator-(time_interp const& ti, double r);
  friend time_interp operator-(double r, time_interp const& ti);

  time_interp & operator*=(double r);
  friend time_interp operator*(time_interp const& ti, double r);
  friend time_interp operator*(double r, time_interp const& ti);

  time_interp & operator/=(double r);
  friend time_interp operator/(time_interp const& ti, double r);
  friend time_interp operator/(double r, time_interp const& ti);

  // Mixed arithmetics time_interp + std::complex<double>
  time_interp & operator+=(std::complex<double> const& z);
  friend time_interp operator+(time_interp const& ti, std::complex<double> const& z);
  friend time_interp operator+(std::complex<double> const& z, time_interp const& ti);

  time_interp & operator-=(std::complex<double> const& z);
  friend time_interp operator-(time_interp const& ti, std::complex<double> const& z);
  friend time_interp operator-(std::complex<double> const& z, time_interp const& ti);

  time_interp & operator*=(std::complex<double> const& z);
  friend time_interp operator*(time_interp const& ti, std::complex<double> const& z);
  friend time_interp operator*(std::complex<double> const& z, time_interp const& ti);

  time_interp & operator/=(std::complex<double> const& z);
  friend time_interp operator/(time_interp const& ti, std::complex<double> const& z);
  friend time_interp operator/(std::complex<double> const& z, time_interp const& ti);

  // Real part
  time_interp real() const {
    return time_interp(mesh, (real_data_type)triqs::arrays::real(data));
  }

  // Imaginary part
  time_interp imag() const {
    return time_interp(mesh, (real_data_type)triqs::arrays::imag(data));
  }

  // Complex conjugate
  time_interp conj() const {
    return time_interp(mesh, (data_type)triqs::arrays::conj(data));
  }
};

}

namespace triqs::utility {

inline bool is_zero(realevol::time_interp const& ti,
                    double tolerance = 100 * std::numeric_limits<double>::epsilon()) {
  return ti.is_zero(tolerance);
}
inline realevol::time_interp conj(realevol::time_interp const& ti) { return ti.conj(); }

}
