# Generated automatically using the command :
# c++2py.py ../c++/solver.hpp -p -m realevol -o realevol --appname realevol --moduledoc "The Real-time evolution solver"
from wrap_generator import *

# The module
module = module_(full_name = "realevol", doc = "The Real-time evolution solver", app_name = "realevol")

# All the triqs C++/Python modules
module.use_module('texpr', 'realevol')
module.use_module('operators', 'realevol')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("solver.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/variant.hpp>
using triqs::operators::many_body_operator;
using namespace realevol;
#include "./converters.hxx"
""")

module.add_enum("ode_solve_method", ["RungeKutta","Lanczos"], "realevol", "Method to solve the Schroedinger equation")

# The class solver
c = class_(
        py_type = "Solver",  # name of the python class
        c_type = "solver",   # name of the C++ class
        doc = r"",   # doc of the C++ class
)

c.add_constructor("""(std::map<std::string,indices_type> gf_struct, std::pair<double,double> time_window, int n_t = 1000)""",
                  doc = """ """)

c.add_method("""void solve (**realevol::solve_parameters_t)""",
             doc = """+----------------+----------------------------+-------------------------------+--------------------------------------------+
| Parameter Name | Type                       | Default                       | Documentation                              |
+================+============================+===============================+============================================+
| h              | operator_t                 |                               | Hamiltonian                                |
+----------------+----------------------------+-------------------------------+--------------------------------------------+
| verbosity      | int                        | 3 on MPI rank 0, 0 otherwise. | Verbosity level                            |
+----------------+----------------------------+-------------------------------+--------------------------------------------+
| hbar           | double                     | 1.0                           | Planck constant                            |
+----------------+----------------------------+-------------------------------+--------------------------------------------+
| method         | realevol::ode_solve_method | Lanczos                       | Method to solve the Schroedinger equation  |
+----------------+----------------------------+-------------------------------+--------------------------------------------+ """)

c.add_property(name = "last_run_parameters",
               getter = cfunction("realevol::solve_parameters_t get_last_run_parameters ()"),
               doc = """Set of parameters used in the last call to solve """)

module.add_class(c)

module.generate_code()
