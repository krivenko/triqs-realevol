from wrap_generator import *

module = module_(full_name = "gf_retime", doc = "Multi-variable Green's functions in real time", app_name = "realevol")

module.add_include("<triqs/gfs.hpp>")
module.use_module('gf', 'triqs')

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_include("<triqs/python_tools/converters/vector.hpp>")
module.add_include("<triqs/python_tools/converters/tuple.hpp>")
module.add_include("<triqs/python_tools/converters/arrays.hpp>")
module.add_include("<triqs/python_tools/converters/h5.hpp>")
module.add_using("namespace triqs::gfs")

def no_ns(s): return s % (("",) * s.count("%s"))
def ns(s): return s % (("triqs::gfs::",) * s.count("%s"))

# Wrapper for gf_mesh<cartesian_product<retime,retime>>
mesh_c_type = "%sgf_mesh<%scartesian_product<%sretime,%sretime>>"
m = class_(py_type = "MeshReTimeReTime",
           c_type = no_ns(mesh_c_type),
           c_type_absolute = ns(mesh_c_type),
           is_printable = True,
           hdf5 = True,
           comparisons = "== !="
          )
m.add_constructor("(gf_mesh<retime> m1, gf_mesh<retime> m2)")
m.add_len(doc = "Size of the mesh")
m.add_iterator(c_cast_type = "std::tuple<double,double>")

m.add_method("triqs::utility::mini_vector<size_t,2> size_of_components()")
m.add_method("std::tuple<gf_mesh<retime>,gf_mesh<retime>> components()")

module.add_class(m)

# Wrapper for g_2t_t
gf_c_type = "%sgf_view<%scartesian_product<%sretime,%sretime>>"
g = class_(py_type = "GfReTime_x_ReTime",
    c_type = no_ns(gf_c_type),
    c_type_absolute = ns(gf_c_type),
    is_printable= True,
    hdf5 = True,
    arithmetic = ("algebra","dcomplex","with_inplace_operators")
   )

g.add_constructor("(%s mesh, mini_vector<size_t,2> shape,"
                  "std::vector<std::vector<std::string>> indices = std::vector<std::vector<std::string>>{},"
                  "std::string name = "")" % no_ns(mesh_c_type),
                  python_precall = "realevol._gf_retime.init")
g.add_member(c_name = "name", c_type = "std::string",  doc = "Name of the Green function (used for plotting only)")
g.add_property(getter = cfunction("%s mesh()" % no_ns(mesh_c_type)), doc ="The mesh"),
g.add_property(name = "data",
               getter = cfunction(calling_pattern="auto result = self_c.data()", signature = "array_view<dcomplex,4>()"),
               doc ="The data ")
g.add_property(name = "target_shape",
               getter = cfunction(calling_pattern="auto result = get_target_shape(self_c)", signature = "shape_type()"),
               doc = "")

module.add_class(g)

module.generate_code()
