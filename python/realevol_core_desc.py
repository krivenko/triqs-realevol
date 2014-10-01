from wrap_generator import *

module = module_(full_name = "pytriqs.applications.realevol.realevol_core", doc = "Time-dependent expression")
module.add_include("<complex>")
module.add_include("<triqs/arrays.hpp>")
module.add_using("namespace realevol")
module.add_using("triqs::utility::is_zero")

# Real time-dependent expression
module.add_include("c++/time_expr_r.hpp")

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
module.add_include("c++/time_expr_c.hpp")

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

## Fermionic operators
module.add_include("<triqs/operators/many_body_operator.hpp>")
module.add_using("namespace triqs::utility")

# The operator class (real-valued)
op_r = class_(
       py_type = "Operator_r",
       c_type = "many_body_operator<time_expr_r>",
       c_type_absolute = "triqs::utility::many_body_operator<realevol::time_expr_r>",
       is_printable = True,
       arithmetic = ("algebra","with_unit","with_unary_minus","time_expr_r","double","std::string")
)

op_r.add_constructor(signature="()", doc="create zero operator")
op_r.add_method("bool is_zero()", doc = "Boolean : is the operator null ?")
module.add_class(op_r)

# Add various overload of c, c_dag to the module Annihilation & Creation operators
for name, doc in [("c","annihilation operator"), ("c_dag","creation operator"), ("n","number operator")]:
    for arg in ("std::string ind", "int ind"):
        module.add_function(name=name+'_r',
                            signature="many_body_operator<time_expr_r> %s<time_expr_r>(%s)"%(name,arg),
                            calling_pattern="auto result = %s<time_expr_r>(ind)"%name,
                            doc=doc)

module.add_function(signature="many_body_operator<time_expr_r> dagger(many_body_operator<time_expr_r> op)",
                    doc="Hermitian conjugate")

# The operator class (complex-valued)
op_c = class_(
       py_type = "Operator_c",
       c_type = "many_body_operator<time_expr_c>",
       c_type_absolute = "triqs::utility::many_body_operator<realevol::time_expr_c>",
       is_printable = True,
       arithmetic = ("algebra","with_unit","with_unary_minus","time_expr_c","double","std::string","std::complex<double>")
)

op_c.add_constructor(signature="()", doc="create zero operator")
op_c.add_method("bool is_zero()", doc = "Boolean : is the operator null ?")
module.add_class(op_c)

# Add various overload of c, c_dag to the module Annihilation & Creation operators
for name, doc in [("c","annihilation operator"), ("c_dag","creation operator"), ("n","number operator")]:
    for arg in ("std::string ind", "int ind"):
        module.add_function(name=name+'_c',
                            signature="many_body_operator<time_expr_c> %s<time_expr_c>(%s)"%(name,arg),
                            calling_pattern="auto result = %s<time_expr_c>(ind)"%name,
                            doc=doc)

module.add_function(signature="many_body_operator<time_expr_c> dagger(many_body_operator<time_expr_c> op)",
                    doc="Hermitian conjugate")

## Solver class
#module.add_include("c++/realevol.hpp")

## Solver class (version for real-valued operators)
#solver_r = class_(
        #py_type = "Solver_r",
        #c_type = "solver<segment_mesh,false>",
        #c_type_absolute = "realevol::solver<segment_mesh,false>",
        #is_printable= False,
        #doc = "Solver class (version for real-valued operators)"
        #)

#solver_r.add_constructor(signature="(std::set<std::string> operator_indices)", doc="create an instance of Solver")
#module.add_class(solver_r)

## Solver class (version for complex-valued operators)
#solver_c = class_(
        #py_type = "Solver_c",
        #c_type = "solver<segment_mesh,true>",
        #c_type_absolute = "realevol::solver<segment_mesh,true>",
        #is_printable= False,
        #doc = "Solver class (version for complex-valued operators)"
        #)

#solver_c.add_constructor(signature="(std::set<std::string> operator_indices)", doc="create an instance of Solver")
#module.add_class(solver_c)

# Generate the module
module.generate_code()