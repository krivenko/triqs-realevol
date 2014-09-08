from wrap_generator import *
from itertools import product

# realevol module
module = module_(full_name = "pytriqs.applications.realevol", doc = "Real-time propagation solver")
module.add_include("<triqs/arrays.hpp>")
module.add_include("c++/time_expr.hpp")
module.add_include("c++/callable_complex.hpp")
module.add_using("namespace realevol")
module.add_using("triqs::utility::is_zero")

# Real time-dependent expression
texpr = class_(
        py_type = "texpr",
        c_type = "time_expr",
        c_type_absolute = "realevol::time_expr",
        is_printable= True,
        arithmetic = ("algebra","with_unit","with_unary_minus","double"),
        doc = "Time-dependent expression (muParser wrapper)"
        )

texpr.add_constructor(signature="(std::string expr)", doc="create expression")
texpr.add_constructor(signature="(double r = 0)", doc="create constant-valued expression")
texpr.add_call(signature="double(double t)", doc="Substitute a time value into the expression")
module.add_class(texpr)

module.add_function(signature="bool is_zero(time_expr te)", doc="Boolean : is zero expression?")
module.add_function(signature="bool is_constant(time_expr te)", doc="Boolean : is constant expression?")

# Complex time-dependent expression
ctexpr = class_(
        py_type = "ctexpr",
        c_type = "callable_complex<time_expr>",
        c_type_absolute = "realevol::callable_complex<realevol::time_expr>",
        is_printable= True,
        arithmetic = ("algebra","with_unit","with_unary_minus","time_expr"),
        doc = "Complex time-dependent expression (muParser wrapper)"
        )

ctexpr_constructor_arg_types = ('double','std::string','time_expr')
for arg_t1, arg_t2 in product(ctexpr_constructor_arg_types,ctexpr_constructor_arg_types):
    ctexpr.add_constructor(signature="(%s re = %s(), %s im = %s())" %((arg_t1,arg_t1,arg_t2,arg_t2)), doc="create expression")
ctexpr.add_constructor(signature="(std::complex<double> z)", doc="create expression")
ctexpr.add_call(signature="std::complex<double>(double t)", doc="Substitute a time value into the expression")
module.add_class(ctexpr)

module.add_function(signature="bool is_zero(callable_complex<time_expr> te)", doc="Boolean : is zero expression?")
module.add_function(signature="bool is_constant(callable_complex<time_expr> te)", doc="Boolean : is constant expression?")

module.generate_code()