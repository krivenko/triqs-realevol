## Time meshes
from wrap_generator import *

module = module_(full_name = "pytriqs.applications.realevol.meshes", doc = "Various time meshes")
module.add_include("<triqs/arrays.hpp>")
module.add_include("c++/any_mesh.hpp")
module.add_using("triqs::gfs::segment_mesh")

# Regular mesh on a segment
rmesh = class_(
        py_type = "rmesh",
        c_type = "segment_mesh",
        c_type_absolute = "triqs::gfs::segment_mesh",
        is_printable= False,
        doc = "Regular time mesh on a segment"
        )

rmesh.add_constructor(signature="(double start, double end, int nodes)", doc="Create regular mesh")
rmesh.add_getitem(signature="double(int node)")
rmesh.add_iterator(c_cast_type="double")
rmesh.add_len(c_name = "size", calling_pattern="int result = self_c.size()")
rmesh.add_property(getter=cfunction("double delta()"), doc="Step of the mesh")
module.add_class(rmesh)

module.generate_code()
