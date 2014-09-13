from wrap_generator import *

module = module_(full_name = "pytriqs.applications.realevol.meshes", doc = "Time meshes")
module.add_include("<triqs/arrays.hpp>")
module.add_include("c++/uniform_mesh.hpp")
module.add_using("namespace realevol")

# Uniform mesh
uniform_mesh = class_(
        py_type = "umesh",
        c_type = "uniform_mesh<std::size_t,double>",
        c_type_absolute = "realevol::uniform_mesh<std::size_t,double>",
        is_printable= False,
        doc = "Unifoirm time mesh"
        )

uniform_mesh.add_constructor(signature="(double start, double end, std::size_t nodes)", doc="create uniform mesh")
uniform_mesh.add_getitem(signature="double(std::size_t node)")
uniform_mesh.add_iterator(c_cast_type="double")
uniform_mesh.add_len(c_name = "size", calling_pattern="int result = self_c.size()")
uniform_mesh.add_property(getter=cfunction("double get_step()"), doc="step of the mesh")
module.add_class(uniform_mesh)

module.generate_code()
