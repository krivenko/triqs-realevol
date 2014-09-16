from wrap_generator import *

# The many_body_operators module
module = module_(full_name = "pytriqs.applications.realevol.operators", doc = "Doc to be written")
module.use_module("texpr")
module.add_include("<triqs/arrays.hpp>")
module.add_include("<triqs/operators/many_body_operator.hpp>")
module.add_using("namespace triqs::utility")
module.add_using("namespace realevol")

# The operator class
op = class_(
        py_type = "Operator",
        c_type = "many_body_operator<time_expr>",
        c_type_absolute = "triqs::utility::many_body_operator<realevol::time_expr>",
        is_printable = True,
        arithmetic = ("algebra","with_unit","with_unary_minus","time_expr","double","std::string")
)

op.add_constructor(signature="()", doc="create zero operator")
op.add_method("bool is_zero()", doc = "Boolean : is the operator null ?")
module.add_class(op)

# Add various overload of c, c_dag to the module Annihilation & Creation operators
for name, doc in [("c","annihilation operator"), ("c_dag","creation operator"), ("n","number operator")]:
    for arg in ("std::string ind", "int ind"):
        module.add_function(name=name,
                            signature="many_body_operator<time_expr> %s<time_expr>(%s)"%(name,arg),
                            calling_pattern="auto result = %s<time_expr>(ind)"%name,
                            doc=doc)

module.generate_code()
