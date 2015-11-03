from wrap_generator import *

module = module_(full_name = "pytriqs.applications.impurity_solvers.realevol.texpr", doc = "Time-dependent expression")
module.add_include("c++/time_expr.hpp")
module.add_using("namespace realevol")
module.add_using("triqs::utility::is_zero")
module.add_using("triqs::utility::conj")

texpr = class_(
        py_type = "texpr",
        c_type = "time_expr",
        c_type_absolute = "realevol::time_expr",
        is_printable = True,
        arithmetic = ("algebra","with_unit","with_unary_minus","double","std::complex<double>"),
        doc = "Time-dependent expression (ExprTk wrapper)"
        )

texpr.add_constructor(signature="()", doc="Create zero expression")
texpr.add_constructor(signature="(double r)", doc="Create real constant-valued expression")
texpr.add_constructor(signature="(std::complex<double> z)", doc="Create complex constant-valued expression")
texpr.add_constructor(signature="(double r, double i)", doc="Create complex constant-valued expression out of two real numbers")
texpr.add_constructor(signature="(std::string re_str)", doc="Create real expression")
texpr.add_constructor(signature="(std::string re_str, std::string im_str)", doc="Create complex expression")
texpr.add_constructor(signature="(double r, std::string im_str)", doc="Create complex expression with a constant real part")
texpr.add_constructor(signature="(std::string re_str, double i)", doc="Create complex expression with a constant imaginary part")

texpr.add_call(signature="std::complex<double>(double t)", doc="Substitute a time value into the expression")
module.add_class(texpr)

module.add_function(signature="bool is_zero(time_expr te)", doc="Boolean : is zero expression?")
module.add_function(signature="bool is_constant(time_expr te)", doc="Boolean : is constant expression?")
module.add_function(signature="time_expr conj(time_expr te)", doc="Complex conjugate of the expression")

module.generate_code()
