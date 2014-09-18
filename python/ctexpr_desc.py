from wrap_generator import *
from itertools import product

# ctexpr module
module = module_(full_name = "pytriqs.applications.realevol.ctexpr", doc = "Complex time-dependent expression")
module.use_module("texpr")
module.add_include("<triqs/arrays.hpp>")
module.add_include("c++/time_expr.hpp")
module.add_include("c++/c_time_expr.hpp")
module.add_using("namespace realevol")
module.add_using("triqs::utility::is_zero")

# Complex time-dependent expression
texpr = class_(
        py_type = "texpr",
        c_type = "c_time_expr",
        c_type_absolute = "realevol::c_time_expr",
        is_printable= True,
        arithmetic = ("algebra","with_unit","with_unary_minus","time_expr","double","std::string"),
        doc = "Complex time-dependent expression (ExprTk wrapper)"
        )

ctexpr_constructor_arg_types = ('double','std::string','time_expr')
for arg_t1, arg_t2 in product(ctexpr_constructor_arg_types,ctexpr_constructor_arg_types):
    texpr.add_constructor(signature="(%s re = %s(), %s im = %s())" %((arg_t1,arg_t1,arg_t2,arg_t2)), doc="create expression")
texpr.add_constructor(signature="(std::complex<double> z)", doc="create expression")
texpr.add_call(signature="std::complex<double>(double t)", doc="Substitute a time value into the expression")
module.add_class(texpr)

module.add_function(signature="bool is_zero(c_time_expr te)", doc="Boolean : is zero expression?")
module.add_function(signature="bool is_constant(c_time_expr te)", doc="Boolean : is constant expression?")

module.generate_code()
