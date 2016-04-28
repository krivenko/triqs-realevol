from wrap_generator import *

module = module_(full_name = "operators", app_name="realevol", doc = "Many-body operator as a function of time")

module.use_module("texpr")

module.add_include("<triqs/operators/many_body_operator.hpp>")
module.add_include("<triqs/python_tools/converters/pair.hpp>")
module.add_include("<triqs/python_tools/converters/vector.hpp>")
module.add_include("<triqs/python_tools/converters/variant_int_string.hpp>")
module.add_include("<triqs/python_tools/converters/h5.hpp>")
module.add_using("namespace triqs::operators")
module.add_using("namespace realevol")

# The operator class
op = class_(
        py_type = "Operator",
        c_type = "many_body_operator_generic<time_expr>",
        c_type_absolute = "triqs::operators::many_body_operator_generic<realevol::time_expr>",
        is_printable= True,
        arithmetic = ("algebra","with_unit","with_unary_minus","realevol::time_expr","double","dcomplex")
        )

op.add_constructor(signature="()", doc="create zero operator")
op.add_constructor(signature="(realevol::time_expr x)", doc="create a constant operator")
op.add_method("bool is_zero()", doc = "Boolean : is the operator null ?")
op.add_iterator(c_cast_type="std::pair<std::vector<std::pair<bool,triqs::operators::indices_t>>, realevol::time_expr>")

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

