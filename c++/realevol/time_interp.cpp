/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2024, I. Krivenko, M. Danilov, P. Kubiczek
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
#include "time_interp.hpp"

namespace realevol {

using namespace std::literals::complex_literals;

//
// Constructors
//

time_interp::time_interp(double r) : interp_t(r) {}
time_interp::time_interp(double r, double i) : interp_t(std::complex<double>(r, i)) {}
time_interp::time_interp(std::complex<double> const& z) : interp_t(z) {}
time_interp::time_interp(mesh_type const& mesh) : interp_t(mesh, .0) {}
time_interp::time_interp(mesh_type const& mesh, double r) : interp_t(mesh, r) {}
time_interp::time_interp(mesh_type const& mesh, double r, double i) : interp_t(mesh, std::complex<double>(r, i)) {}
time_interp::time_interp(mesh_type const& mesh, std::complex<double> const& z) : interp_t(mesh, z) {}
time_interp::time_interp(mesh_type const& mesh, real_data_type const& r) : interp_t(mesh, r) {}
time_interp::time_interp(mesh_type const& mesh, real_data_type const& r, double i) : interp_t(mesh, r) { data += 1i*i; }
time_interp::time_interp(mesh_type const& mesh, double r, real_data_type const& i) : interp_t(mesh, 1i*i) { data += r; }
time_interp::time_interp(mesh_type const& mesh, real_data_type const& r, real_data_type const& i) : interp_t(mesh, r + 1i*i) {}
time_interp::time_interp(mesh_type const& mesh, data_type const& z) : interp_t(mesh, z) {}

//
// Arithmetics
//

time_interp time_interp::operator-() const {
  return is_constant_ ? time_interp(mesh, -data(0)) : time_interp(mesh, -data);
}

time_interp & time_interp::operator+=(const time_interp& ti) {
  if(is_constant_) {
    if(ti.is_constant_)
      data(0) += ti.data(0);
    else {
      auto z = data(0);
      *this = ti;
      data += z;
    }
  } else {
    if(ti.is_constant_)
      data += ti.data(0);
    else {
      TRIQS_ASSERT(mesh == ti.mesh);
      data += ti.data;
    }
  }
  return *this;
}

time_interp & time_interp::operator-=(const time_interp& ti) {
  if(is_constant_) {
    if(ti.is_constant_)
      data(0) -= ti.data(0);
    else {
      auto z = data(0);
      *this = ti;
      data -= z;
    }
  } else {
    if(ti.is_constant_)
      data -= ti.data(0);
    else {
      TRIQS_ASSERT(mesh == ti.mesh);
      data -= ti.data;
    }
  }
  return *this;
}

time_interp & time_interp::operator*=(const time_interp& ti) {
  if(is_constant_) {
    if(ti.is_constant_)
      data(0) *= ti.data(0);
    else {
      auto z = data(0);
      *this = ti;
      data *= z;
    }
  } else {
    if(ti.is_constant_)
      data *= ti.data(0);
    else {
      TRIQS_ASSERT(mesh == ti.mesh);
      data *= ti.data;
    }
  }
  return *this;
}

time_interp & time_interp::operator/=(const time_interp& ti) {
  if(is_constant_) {
    if(ti.is_constant_)
      data(0) /= ti.data(0);
    else {
      auto z = data(0);
      *this = ti;
      data /= z;
    }
  } else {
    if(ti.is_constant_)
      data /= ti.data(0);
    else {
      TRIQS_ASSERT(mesh == ti.mesh);
      data /= ti.data;
    }
  }
  return *this;
}

//
// Mixed arithmetics time_interp + double
//

time_interp & time_interp::operator+=(double r) {
  if(is_constant_)
    data(0) += r;
  else
    data() += r;
  return *this;
}

time_interp operator+(time_interp const& ti, double r) {
  time_interp res(ti); res += r;
  return res;
}

time_interp operator+(double r, time_interp const& ti) {
  time_interp res(ti); res += r;
  return res;
}

time_interp & time_interp::operator-=(double r) {
  if(is_constant_)
    data(0) -= r;
  else
    data() -= r;
  return *this;
}

time_interp operator-(time_interp const& ti, double r) {
  time_interp res(ti); res -= r;
  return res;
}

time_interp operator-(double r, time_interp const& ti) {
  time_interp res(-ti); res += r;
  return res;
}

time_interp & time_interp::operator*=(double r) {
  if(is_constant_)
    data(0) *= r;
  else
    data() *= r;
  return *this;
}

time_interp operator*(time_interp const& ti, double r) {
  time_interp res(ti); res *= r;
  return res;
}

time_interp operator*(double r, time_interp const& ti) {
  time_interp res(ti); res *= r;
  return res;
}

time_interp & time_interp::operator/=(double r) {
  if(is_constant_)
    data(0) /= r;
  else
    data() /= r;
  return *this;
}

time_interp operator/(time_interp const& ti, double r) {
  time_interp res(ti); res /= r;
  return res;
}

time_interp operator/(double r, time_interp const& ti) {
  return time_interp(ti.mesh, r / ti.data);
}

//
// Mixed arithmetics time_interp + std::complex<double>
//

time_interp & time_interp::operator+=(std::complex<double> const& z) {
  if(is_constant_)
    data(0) += z;
  else
    data() += z;
  return *this;
}

time_interp operator+(time_interp const& ti, std::complex<double> const& z) {
  time_interp res(ti); res += z;
  return res;
}

time_interp operator+(std::complex<double> const& z, time_interp const& ti) {
  time_interp res(ti); res += z;
  return res;
}

time_interp & time_interp::operator-=(std::complex<double> const& z) {
  if(is_constant_)
    data(0) -= z;
  else
    data() -= z;
  return *this;
}

time_interp operator-(time_interp const& ti, std::complex<double> const& z) {
  time_interp res(ti); res -= z;
  return res;
}

time_interp operator-(std::complex<double> const& z, time_interp const& ti) {
  time_interp res(-ti); res += z;
  return res;
}

time_interp & time_interp::operator*=(std::complex<double> const& z) {
  if(is_constant_)
    data(0) *= z;
  else
    data() *= z;
  return *this;
}

time_interp operator*(time_interp const& ti, std::complex<double> const& z) {
  time_interp res(ti); res *= z;
  return res;
}

time_interp operator*(std::complex<double> const& z, time_interp const& ti) {
  time_interp res(ti); res *= z;
  return res;
}

time_interp & time_interp::operator/=(std::complex<double> const& z) {
  if(is_constant_)
    data(0) /= z;
  else
    data() /= z;
  return *this;
}

time_interp operator/(time_interp const& ti, std::complex<double> const& z) {
  time_interp res(ti); res /= z;
  return res;
}

time_interp operator/(std::complex<double> const& z, time_interp const& ti) {
  return time_interp(ti.mesh, z / ti.data);
}

}
