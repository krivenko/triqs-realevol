from wrap_generator import *

# texpr module
module = module_(full_name = "pytriqs.applications.realevol.realevol_core", doc = "Time-dependent expression")
module.add_include("<complex>")
module.add_include("<triqs/arrays.hpp>")
module.add_include("c++/time_expr_r.hpp")
module.add_include("c++/time_expr_c.hpp")
module.add_using("namespace realevol")
module.add_using("triqs::utility::is_zero")

# Real time-dependent expression
texpr_r = class_(
        py_type = "texpr_r",
        c_type = "time_expr_r",
        c_type_absolute = "realevol::time_expr_r",
        is_printable= True,
        arithmetic = ("algebra","with_unit","with_unary_minus","double","std::string"),
        doc = "Real-valued time-dependent expression (ExprTk wrapper)"
        )

texpr_r.add_constructor(signature="(std::string expr)", doc="create expression")
texpr_r.add_constructor(signature="(double r = 0)", doc="create constant-valued expression")
texpr_r.add_call(signature="double(double t)", doc="Substitute a time value into the expression")
module.add_class(texpr_r)

module.add_function(signature="bool is_zero(time_expr_r te)", doc="Boolean : is zero expression?")
module.add_function(signature="bool is_constant(time_expr_r te)", doc="Boolean : is constant expression?")

# Complex time-dependent expression
texpr_c = class_(
        py_type = "texpr_c",
        c_type = "time_expr_c",
        c_type_absolute = "realevol::time_expr_c",
        is_printable= True,
        arithmetic = ("algebra","with_unit","with_unary_minus","time_expr_r","double","std::string","std::complex<double>"),
        doc = "Complex-valued time-dependent expression (ExprTk wrapper)"
        )

from itertools import product

constructor_arg_types = ('double','std::string','time_expr_r')
for arg_t1, arg_t2 in product(constructor_arg_types,constructor_arg_types):
    texpr_c.add_constructor(signature="(%s re = %s(), %s im = %s())" %((arg_t1,arg_t1,arg_t2,arg_t2)), doc="create expression")
texpr_c.add_constructor(signature="(std::complex<double> z)", doc="create expression")
texpr_c.add_call(signature="std::complex<double>(double t)", doc="Substitute a time value into the expression")
module.add_class(texpr_c)

module.add_function(signature="bool is_zero(time_expr_c te)", doc="Boolean : is zero expression?")
module.add_function(signature="bool is_constant(time_expr_c te)", doc="Boolean : is constant expression?")

## Time meshes
module.add_include("<triqs/gfs/meshes/segment.hpp>")
module.add_using("triqs::gfs::segment_mesh")

# Regular mesh on a segment
rmesh = class_(
        py_type = "rmesh",
        c_type = "segment_mesh",
        c_type_absolute = "triqs::gfs::segment_mesh",
        is_printable= False,
        doc = "Regular time mesh on a segment"
        )

rmesh.add_constructor(signature="(double start, double end, int nodes)", doc="create uniform mesh")
rmesh.add_getitem(signature="double(int node)")
rmesh.add_iterator(c_cast_type="double")
rmesh.add_len(c_name = "size", calling_pattern="int result = self_c.size()")
rmesh.add_property(getter=cfunction("double delta()"), doc="step of the mesh")
module.add_class(rmesh)

# Generate the module
module.generate_code()