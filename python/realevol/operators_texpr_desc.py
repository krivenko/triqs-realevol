# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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

module = module_(full_name = "operators_texpr", app_name = "realevol", doc = "Many-body operator with TExpr coefficients")

module.add_imports('realevol.texpr')

module.add_include("<triqs/utility/variant_extensions.hpp>")
module.add_include("<realevol/time_expr.hpp>")
module.add_include("<realevol/operators/many_body_operator.hpp>")

# FIXME: This include makes a little sense here, but without it the module will
# SEGFAULT in `cpp2py::py_converter<std::complex<double>>::py2c()` whenever
# a complex value is passed from Python to a C++ function.
module.add_include("<triqs/arrays.hpp>")
module.add_include("<triqs/cpp2py_converters.hpp>")

module.add_preamble("""
using namespace realevol::operators;
using namespace realevol;
""")

# The operator class
op = class_(
        py_type = "Operator",
        c_type = "many_body_operator_generic<time_expr>",
        c_type_absolute = "realevol::operators::many_body_operator_generic<realevol::time_expr>",
        is_printable= True,
        arithmetic = ("algebra","with_unit","with_unary_minus","realevol::time_expr","double","dcomplex")
        )

op.add_constructor(signature="()", doc="create zero operator")
op.add_constructor(signature="(realevol::time_expr x)", doc="create a constant operator")
op.add_method("bool is_zero()", doc = "Boolean : is the operator null ?")
op.add_iterator(c_cast_type="std::pair<std::vector<std::pair<bool,indices_t>>, realevol::time_expr>")

module.add_class(op)

# Annihilation & Creation operators
for name, doc in [("c","Fermionic annihilation operator"),
                  ("c_dag","Fermionic creation operator"),
                  ("n","Fermionic number operator"),
                  ("a","Bosonic annihilation operator"),
                  ("a_dag","Bosonic creation operator")] :
    for args in [[("","")],
            [("std::string","ind1")],
            [("std::string","ind1"),("std::string","ind2")],
            [("int","i"),("std::string","ind1")],
            [("std::string","ind1"),("int","i")],
            [("int","i"), ("int","j")]
            ]:

        signature = "many_body_operator_generic<time_expr>(" +','.join(map(lambda a:"%s %s"%(a[0],a[1]), args)) + ")"
        calling_pattern = "auto result = %s<time_expr>("%name +','.join(map(lambda a: str(a[1]), args)) + ")"
        module.add_function(name = name, signature = signature, calling_pattern = calling_pattern, doc = doc)

module.add_function("many_body_operator_generic<time_expr> dagger(many_body_operator_generic<time_expr> Op)",
                    doc = "Return the Hermitian conjugate of the operator")

module.generate_code()
