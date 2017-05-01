# Generated automatically using the command :
# c++2py.py ../c++/solver.hpp -p -m realevol -o realevol --appname realevol --moduledoc "The Real-time evolution solver"
from wrap_generator import *

# The module
module = module_(full_name = "realevol", doc = "The Real-time evolution solver", app_name = "realevol")

# All the triqs C++/Python modules
module.use_module('texpr', 'realevol')
module.use_module('operators', 'realevol')
module.use_module('init_state', 'realevol')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("solver.hpp")
module.add_include("utility.hpp")
module.add_include("<triqs/gfs/gf_tests.hpp>")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_include("<triqs/python_tools/converters/pair.hpp>")
module.add_include("<triqs/python_tools/converters/map.hpp>")
module.add_include("<triqs/python_tools/converters/vector.hpp>")
module.add_include("<triqs/python_tools/converters/gf.hpp>")
module.add_include("<triqs/python_tools/converters/block_gf.hpp>")
module.add_using("realevol::operators::many_body_operator") # FIXME
module.add_using("namespace triqs::gfs")
module.add_preamble("""
using namespace realevol;
#include "init_state_converters.hxx"
#include "compute_2t_obs_converters.hxx"
""")

module.add_enum("h_interpolation", ["Rectangle", "Trapezoid", "Simpson"], "realevol",
                "Hamiltonian interpolation between time slices")

# The class solver
c = class_(
        py_type = "Solver",  # name of the python class
        c_type = "solver",   # name of the C++ class
        c_type_absolute = "realevol::solver",
        doc = "The Real-time evolution solver",   # doc of the C++ class
)

c.add_constructor("""(gf_struct_t gf_struct, chi_indices_t chi_indices, double t_max, int n_t = 1000)""",
                  doc = """ """)

c.add_method("""void compute_2t_obs (**compute_2t_obs_parameters_t)""",
             doc = """+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+
| Parameter Name          | Type                      | Default                       | Documentation                                                                 |
+=========================+===========================+===============================+===============================================================================+
| h                       | operator_t                | --                            | Hamiltonian                                                                   |
+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+
| verbosity               | int                       | 3 on MPI rank 0, 0 otherwise. | Verbosity level                                                               |
+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+
| hbar                    | double                    | 1.0                           | Planck constant                                                               |
+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+
| hamiltonian_interpol    | realevol::h_interpolation | Rectangle                     | Hamiltonian interpolation between time slices                                 |
+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+
| lanczos_min_matrix_size | int                       | 11                            | Use Lanczos algorithm to exponentiate matrices of this size or bigger         |
+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+
| lanczos_gs_energy_tol   | std::map<long, double>    | {}                            | Lanczos convergence threshold for the GS energy, for each invariant subspace  |
+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+
| lanczos_max_krylov_dim  | std::map<long, int>       | {}                            | Maximal dimension of the Krylov space, for each invariant subspace            |
+-------------------------+---------------------------+-------------------------------+-------------------------------------------------------------------------------+""")

c.add_property(name = "initial_state",
               getter = cfunction("init_state const& get_initial_state()"),
               setter = cfunction("void set_initial_state(init_state initial_state)"),
               doc = """Initial state at t=0""")

c.add_property(name = "last_compute_2t_obs_parameters",
               getter = cfunction("compute_2t_obs_parameters_t get_last_compute_2t_obs_parameters ()"),
               doc = """Set of parameters used in the last call to solve""")

c.add_property(name = "g_l",
               getter = cfunction("block_gf_2t_view get_g_l ()"),
               doc = """Lesser GF in real time""")

c.add_property(name = "g_g",
               getter = cfunction("block_gf_2t_view get_g_g ()"),
               doc = """Greater GF in real time""")

c.add_property(name = "chi",
               getter = cfunction("gf_2t_view get_chi ()"),
               doc = """Susceptibility in real time""")

module.add_class(c)

module.add_function("std::pair<block_gf_2t_view,block_gf_2t_view> make_gf_ret_adv(block_gf_2t_view g_l, block_gf_2t_view g_g)",
                    doc = """Compute retarded and advanced Green's functions out of the lesser and greater components""")

# Comparison tests
module.add_function(name = "assert_gfs_are_close",
                    signature = "void(gf_2t_view x, gf_2t_view y, double precision=1.e-6)",
                    doc = """Compare two real time GFs""")

module.add_function(name = "assert_block_gfs_are_close",
                    signature = "void(block_gf_2t_view x, block_gf_2t_view y, double precision=1.e-6)",
                    doc = """Compare two real time block GFs""")

module.generate_code()
