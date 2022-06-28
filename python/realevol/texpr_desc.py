# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2022, I. Krivenko, M. Danilov, P. Kubiczek
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

module = module_(full_name = "texpr", app_name = "realevol", doc = "Time-dependent expression")

module.add_include("<realevol/time_expr.hpp>")

# FIXME: This include makes a little sense here, but without it the module will
# SEGFAULT in `cpp2py::py_converter<std::complex<double>>::py2c()` whenever
# a complex value is passed from Python to a C++ function.
module.add_include("<triqs/arrays.hpp>")
module.add_include("<triqs/cpp2py_converters.hpp>")

module.add_preamble("""
using namespace realevol;
using triqs::utility::is_zero;
using triqs::utility::conj;
""")

c = class_(
    py_type = "TExpr",
    c_type = "time_expr",
    c_type_absolute = "realevol::time_expr",
    is_printable = True,
    arithmetic = ("algebra","with_unit","with_unary_minus","double","std::complex<double>"),
    doc = "Time-dependent expression (ExprTk wrapper)"
    )

c.add_constructor(signature="()", doc="Create zero expression")
c.add_constructor(signature="(double r)", doc="Create real constant-valued expression")
c.add_constructor(signature="(std::complex<double> z)", doc="Create complex constant-valued expression")
c.add_constructor(signature="(double r, double i)", doc="Create complex constant-valued expression out of two real numbers")
c.add_constructor(signature="(std::string re_str)", doc="Create real expression")
c.add_constructor(signature="(std::string re_str, std::string im_str)", doc="Create complex expression")
c.add_constructor(signature="(double r, std::string im_str)", doc="Create complex expression with a constant real part")
c.add_constructor(signature="(std::string re_str, double i)", doc="Create complex expression with a constant imaginary part")

c.add_call(signature="std::complex<double>(double t)", doc="Substitute a time value into the expression")

module.add_class(c)

module.add_function(signature="bool is_zero(time_expr te)", doc="Boolean : is zero expression?")
module.add_function(signature="bool is_constant(time_expr te)", doc="Boolean : is constant expression?")
module.add_function(signature="time_expr conj(time_expr te)", doc="Complex conjugate of the expression")

module.generate_code()
