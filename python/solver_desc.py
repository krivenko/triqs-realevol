from wrap_generator import *

module = module_(full_name = "pytriqs.applications.realevol.solver", doc = "Real-time evolution solver")
module.add_include("<triqs/arrays.hpp>")
module.add_include("c++/uniform_mesh.hpp")
module.add_include("c++/realevol.hpp")
module.add_using("namespace realevol")

# Solver class (version for real-valued operators)
solver = class_(
        py_type = "Solver",
        c_type = "solver<uniform_mesh<>,false>",
        c_type_absolute = "realevol::solver<realevol::uniform_mesh<>,false>",
        is_printable= False,
        doc = "Solver class (version for real-valued operators)"
        )

solver.add_constructor(signature="(std::vector<std::string> operator_indices)", doc="create an instance of Solver")
module.add_class(solver)

# Solver class (version for complex-valued operators)
csolver = class_(
        py_type = "SolverComplex",
        c_type = "solver<uniform_mesh<>,true>",
        c_type_absolute = "realevol::solver<realevol::uniform_mesh<>,true>",
        is_printable= False,
        doc = "Solver class (version for complex-valued operators)"
        )

csolver.add_constructor(signature="(std::vector<std::string> operator_indices)", doc="create an instance of Solver")
module.add_class(csolver)

module.generate_code()
