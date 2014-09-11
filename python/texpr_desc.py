from wrap_generator import *

# realevol module
module = module_(full_name = "pytriqs.applications.realevol.texpr", doc = "Time-dependent expression")
module.add_include("<triqs/arrays.hpp>")
module.add_include("c++/time_expr.hpp")
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

module.generate_code()