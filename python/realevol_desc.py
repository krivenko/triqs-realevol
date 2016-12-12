# Generated automatically using the command :
# c++2py.py ../c++/solver.hpp -p -m realevol -o realevol --appname realevol --moduledoc "The Real-time evolution solver"
from wrap_generator import *

# The module
module = module_(full_name = "realevol", doc = "The Real-time evolution solver", app_name = "realevol")

# All the triqs C++/Python modules
module.use_module('texpr', 'realevol')
module.use_module('operators', 'realevol')
module.use_module('gf_retime', 'realevol')
module.use_module('init_state', 'realevol')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("solver.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_include("<triqs/python_tools/converters/pair.hpp>")
module.add_include("<triqs/python_tools/converters/map.hpp>")
module.add_include("<triqs/python_tools/converters/vector.hpp>")
module.add_using("realevol::operators::many_body_operator") # FIXME
module.add_using("namespace triqs::gfs")
module.add_preamble("""
using namespace realevol;
#include "init_state_converters.hxx"
#include "compute_gf_converters.hxx"
""")

# The class solver
c = class_(
        py_type = "Solver",  # name of the python class
        c_type = "solver",   # name of the C++ class
        c_type_absolute = "realevol::solver",
        doc = "The Real-time evolution solver",   # doc of the C++ class
)

c.add_constructor("""(std::map<std::string,indices_type> gf_struct, std::pair<double,double> time_window, int n_t = 1000)""",
                  doc = """ """)

c.add_method("""void compute_gf (**compute_gf_parameters_t)""",
             doc = """+----------------+----------------------+-------------------------------+-----------------------+
| Parameter Name | Type                 | Default                       | Documentation         |
+================+======================+===============================+=======================+
| h              | operator_t           | --                            | Hamiltonian           |
+----------------+----------------------+-------------------------------+-----------------------+
| initial_state  | realevol::init_state | --                            | Initial state at t=0  |
+----------------+----------------------+-------------------------------+-----------------------+
| verbosity      | int                  | 3 on MPI rank 0, 0 otherwise. | Verbosity level       |
+----------------+----------------------+-------------------------------+-----------------------+
| hbar           | double               | 1.0                           | Planck constant       |
+----------------+----------------------+-------------------------------+-----------------------+""")

c.add_property(name = "initial_state",
               getter = cfunction("init_state const& get_initial_state()"),
               setter = cfunction("void set_initial_state(init_state initial_state)"),
               doc = """Initial state at t=0""")

c.add_property(name = "last_compute_gf_parameters",
               getter = cfunction("compute_gf_parameters_t get_last_compute_gf_parameters ()"),
               doc = """Set of parameters used in the last call to solve""")

c.add_property(name = "g_l",
               getter = cfunction("block_gf_view<cartesian_product<retime, retime>> get_g_l ()"),
               doc = """Lesser GF in real time""")

c.add_property(name = "g_g",
               getter = cfunction("block_gf_view<cartesian_product<retime, retime>> get_g_g ()"),
               doc = """Greater GF in real time""")

c.add_property(name = "g_ret",
               getter = cfunction("block_gf_view<cartesian_product<retime, retime>> get_g_ret ()"),
               doc = """Retarded GF in real time""")

c.add_property(name = "g_adv",
               getter = cfunction("block_gf_view<cartesian_product<retime, retime>> get_g_adv ()"),
               doc = """Advanced GF in real time""")

module.add_class(c)

module.generate_code()
