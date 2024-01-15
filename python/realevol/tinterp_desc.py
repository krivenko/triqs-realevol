# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2024, I. Krivenko, M. Danilov, P. Kubiczek
#
# realevol is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# realevol is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# realevol (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
# ##############################################################################

from cpp2py.wrap_generator import *

module = module_(full_name = "tinterp", app_name = "realevol", doc = "Linear interpolator on time mesh")

module.add_imports('triqs.gf')

module.add_include("<realevol/time_interp.hpp>")

module.add_include("<triqs/cpp2py_converters.hpp>")

module.add_preamble("""
using namespace realevol;
using namespace triqs;
using triqs::utility::is_zero;
using triqs::utility::conj;
""")

c = class_(
    py_type = "TInterp",
    c_type = "time_interp",
    c_type_absolute = "realevol::time_interp",
    is_printable = True,
    arithmetic = ("algebra","with_unit","with_unary_minus","double","std::complex<double>"),
    comparisons = "==",
    doc = "Linear interpolator on a time mesh"
    )

c.add_constructor(signature="()", doc="Create zero interpolator")
c.add_constructor(signature="(mesh::retime m)", doc="Create zero interpolator")
c.add_constructor(signature="(double r)", doc="Create real constant interpolator")
c.add_constructor(signature="(mesh::retime m, double r)", doc="Create real constant interpolator")
c.add_constructor(signature="(std::complex<double> z)", doc="Create complex constant interpolator")
c.add_constructor(signature="(mesh::retime m, std::complex<double> z)", doc="Create complex constant interpolator")
c.add_constructor(signature="(double r, double i)", doc="Create complex constant interpolator out of two real numbers")
c.add_constructor(signature="(mesh::retime m, double r, double i)", doc="Create complex constant interpolator out of two real numbers")
c.add_constructor(signature="(mesh::retime m, nda::array<double, 1> r)", doc="Create real interpolator")
c.add_constructor(signature="(mesh::retime m, nda::array<double, 1> r, double i)",
                  doc="Create complex interpolator with a constant imaginary part")
c.add_constructor(signature="(mesh::retime m, double r, nda::array<double, 1> i)",
                  doc="Create complex interpolator with a constant real part")
c.add_constructor(signature="(mesh::retime m, nda::array<double, 1> r, nda::array<double, 1> i)",
                  doc="Create complex interpolator")
c.add_constructor(signature="(mesh::retime m, nda::array<std::complex<double>, 1> z)",
                  doc="Create complex interpolator")

c.add_call(signature="std::complex<double>(double t)", doc="Interpolation result at time point t")
c.add_property(getter=cfunction("nda::array<std::complex<double>, 1> get_data()"),
               name="data",
               doc="Values at interpolation nodes")
module.add_class(c)

module.add_function(signature="bool is_zero(time_interp ti)", doc="Boolean : is zero interpolator?")
module.add_function(signature="bool is_constant(time_interp ti)", doc="Boolean : is constant interpolator?")
module.add_function(signature="time_interp conj(time_interp ti)", doc="Complex conjugate of the interpolator")

module.generate_code()
