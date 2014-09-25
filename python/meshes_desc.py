from wrap_generator import *

module = module_(full_name = "pytriqs.applications.realevol.meshes", doc = "Time meshes")
module.add_include("<triqs/arrays.hpp>")
module.add_include("<triqs/gfs/meshes/segment.hpp>")
module.add_using("triqs::gfs::segment_mesh")

# Uniform mesh
uniform_mesh = class_(
        py_type = "umesh",
        c_type = "segment_mesh",
        c_type_absolute = "triqs::gfs::segment_mesh",
        is_printable= False,
        doc = "Linear time mesh on a segment"
        )

uniform_mesh.add_constructor(signature="(double start, double end, int nodes)", doc="create uniform mesh")
uniform_mesh.add_getitem(signature="double(int node)")
uniform_mesh.add_iterator(c_cast_type="double")
uniform_mesh.add_len(c_name = "size", calling_pattern="int result = self_c.size()")
uniform_mesh.add_property(getter=cfunction("double delta()"), doc="step of the mesh")
module.add_class(uniform_mesh)

module.generate_code()
