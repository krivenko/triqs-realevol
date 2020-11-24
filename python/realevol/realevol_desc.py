# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
#
# realevol is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# realevol is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# realevol (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
# ##############################################################################

from cpp2py.wrap_generator import *

module = module_(full_name = "realevol", doc = "The Real-time evolution solver", app_name = "realevol")

module.add_imports('triqs.gf',
                   'realevol.texpr', 'realevol.tinterp',
                   'realevol.operators_texpr', 'realevol.operators_tinterp',
                   'realevol.init_state')

module.add_include("<realevol/solver.hpp>")
module.add_include("<realevol/utility.hpp>")
module.add_include("<triqs/gfs/gf_tests.hpp>")

module.add_include("<triqs/cpp2py_converters.hpp>")
module.add_include("<cpp2py/converters/pair.hpp>")
module.add_include("<cpp2py/converters/map.hpp>")
module.add_include("<cpp2py/converters/vector.hpp>")
module.add_include("<cpp2py/converters/variant.hpp>")

module.add_preamble("""
using namespace triqs::gfs;
using namespace realevol;
using realevol::operators::many_body_operator;
""")

module.add_enum("h_interpolation", ["Rectangle", "Trapezoid", "Simpson"], "realevol",
                "Hamiltonian interpolation between time slices")

# The class solver
c = class_(
    py_type = "Solver",  # name of the python class
    c_type = "solver",   # name of the C++ class
    c_type_absolute = "realevol::solver",
    doc = "The Real-time evolution solver"   # doc of the C++ class
)

c.add_constructor("""(gf_struct_t gf_struct, chi_indices_t chi_indices, double t_max, int n_t = 1000)""",
                  doc = """ """)

compute_2t_obs_doc = """\
h (%s) - time-dependent Hamiltonian
params (dict) - solver parameters

+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| Parameter Name          | Type                              | Default                      | Documentation                                                                        |
+=========================+===================================+==============================+======================================================================================+
| verbosity               | int                               | 3 on MPI rank 0, 0 otherwise | Verbosity level                                                                      |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| compute_g_l             | bool                              | True                         | Compute lesser Green's function                                                      |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| compute_g_g             | bool                              | True                         | Compute greater Green's function                                                     |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| compute_chi             | bool                              | True                         | Compute susceptibility                                                               |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| t_range                 | (float, float)                    | (-inf, inf)                  | Compute components of observables with the first time argument within this range     |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| tp_range                | (float, float)                    | (-inf, inf)                  | Compute components of observables with the second time argument within this range    |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| delta_t_max             | float                             | inf                          | Compute components of observables with time arguments satisfying |t-t'|<=delta_t_max |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| hbar                    | float                             | 1.0                          | Planck constant                                                                      |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| hamiltonian_interpol    | str: Rectangle|Trapezoid|Simpson  | Rectangle                    | Hamiltonian interpolation between time slices                                        |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| lanczos_min_matrix_size | int                               | 11                           | Use Lanczos algorithm to exponentiate matrices of this size or bigger                |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| lanczos_gs_energy_tol   | dict(int : float)                 | {}                           | Lanczos convergence threshold for the GS energy, for each invariant subspace         |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
| lanczos_max_krylov_dim  | dict(int : int)                   | {}                           | Maximal dimension of the Krylov space, for each invariant subspace                   |
+-------------------------+-----------------------------------+------------------------------+--------------------------------------------------------------------------------------+
"""

c.add_method("""void compute_2t_obs (time_expr_operator_t h, compute_2t_obs_parameters_t params)""",
             doc = compute_2t_obs_doc % "operators_texpr.Operator")
c.add_method("""void compute_2t_obs (time_interp_operator_t h, compute_2t_obs_parameters_t params)""",
             doc = compute_2t_obs_doc % "operators_tinterp.Operator")

conv = converter_(
    c_type = "realevol::compute_2t_obs_parameters_t",
    doc = r"""Parameters of compute_2t_obs()""",
)

conv.add_member(c_name = "verbosity",
                c_type = "int",
                initializer = "((mpi::communicator().rank() == 0) ? 3 : 0)",
                doc = r"""Verbosity level""")

conv.add_member(c_name = "compute_g_l",
                c_type = "bool",
                initializer = "true",
                doc = r"""Compute lesser Green's function""")

conv.add_member(c_name = "compute_g_g",
                c_type = "bool",
                initializer = "true",
                doc = r"""Compute greater Green's function""")

conv.add_member(c_name = "compute_chi",
                c_type = "bool",
                initializer = "true",
                doc = r"""Compute susceptibility""")

conv.add_member(c_name = "t_range",
                c_type = "std::pair<double, double>",
                initializer = "std::pair<double, double>{-INFINITY, INFINITY}",
                doc = r"""Compute components of observables with the first time argument within this range""")

conv.add_member(c_name = "tp_range",
                c_type = "std::pair<double, double>",
                initializer = "std::pair<double, double>{-INFINITY, INFINITY}",
                doc = r"""Compute components of observables with the second time argument within this range""")

conv.add_member(c_name = "delta_t_max",
                c_type = "double",
                initializer = "INFINITY",
                doc = r"""Compute components of observables with time arguments satisfying |t-t'|<=delta_t_max""")

conv.add_member(c_name = "hbar",
                c_type = "double",
                initializer = "1.0",
                doc = r"""Planck's constant""")

conv.add_member(c_name = "hamiltonian_interpol",
                c_type = "h_interpolation",
                initializer = "h_interpolation::Rectangle",
                doc = r"""Hamiltonian interpolation between time slices""")

conv.add_member(c_name = "lanczos_min_matrix_size",
                c_type = "int",
                initializer = "11",
                doc = r"""Use Lanczos algorithm to exponentiate matrices of this size or bigger""")

conv.add_member(c_name = "lanczos_gs_energy_tol",
                c_type = "std::map<long, double>",
                initializer = "{}",
                doc = r"""Lanczos convergence threshold for the GS energy, for each invariant subspace""")

conv.add_member(c_name = "lanczos_max_krylov_dim",
                c_type = "std::map<long, int>",
                initializer = "{}",
                doc = r"""Maximal dimension of the Krylov space, for each invariant subspace""")

module.add_converter(conv)

c.add_method("void set_initial_state(init_state initial_state)", doc = """Set initial state at t=0""")

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

## Comparison tests
module.add_function(name = "assert_gfs_are_close",
                    signature = "void(gf_2t_view x, gf_2t_view y, double precision=1.e-6)",
                    doc = """Compare two real time GFs""")

module.add_function(name = "assert_block_gfs_are_close",
                    signature = "void(block_gf_2t_view x, block_gf_2t_view y, double precision=1.e-6)",
                    doc = """Compare two real time block GFs""")

module.generate_code()
