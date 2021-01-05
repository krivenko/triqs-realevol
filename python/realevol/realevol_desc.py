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

module.add_include("<realevol/compute.hpp>")
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
using time_expr_operator_t = realevol::operators::many_body_operator_generic<time_expr>;
using time_interp_operator_t = realevol::operators::many_body_operator_generic<time_interp>;
""")

module.add_enum("h_interpolation", ["Rectangle", "Trapezoid", "Simpson"], "realevol",
                "Hamiltonian interpolation between time slices")

#
# compute_expectval()
#

module.add_function(
    name = "compute_expectval",
    signature = "expectval_container_t(static_operator_t op,"
                "init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh,"
                "solver_parameters_t params)",
    calling_pattern = "auto result = compute_expectval(op, initial_state, h, t_mesh, params)",
    doc = """Compute expectation value of operator 'op' as a function of time"""
)

module.add_function(
    name = "compute_expectval",
    signature = "expectval_container_t(static_operator_t op,"
                "init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh,"
                "solver_parameters_t params)",
    calling_pattern = "auto result = compute_expectval(op, initial_state, h, t_mesh, params)",
    doc = """Compute expectation value of operator 'op' as a function of time"""
)

#
# compute_correlator_2t()
#

module.add_function(
    name = "compute_correlator_2t",
    signature = "correlator_2t_container_t(static_operator_t op1, static_operator_t op2,"
                "init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh,"
                "solver_parameters_t params)",
    calling_pattern = "auto result = compute_correlator_2t(op1, op2, initial_state, h, t_mesh, params)",
    doc = """Compute a 2-point correlator of operators 'op1' and 'op2' with their time arguments defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_correlator_2t",
    signature = "correlator_2t_container_t(static_operator_t op1, static_operator_t op2,"
                "init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh,"
                "solver_parameters_t params)",
    calling_pattern = "auto result = compute_correlator_2t(op1, op2, initial_state, h, t_mesh, params)",
    doc = """Compute a 2-point correlator of operators 'op1' and 'op2' with their time arguments defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

#
# compute_correlator_3t()
#

module.add_function(
    name = "compute_correlator_3t",
    signature = "correlator_3t_container_t(static_operator_t op1, static_operator_t op2, static_operator_t op3,"
                "init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh,"
                "solver_parameters_t params)",
    calling_pattern = "auto result = compute_correlator_3t(op1, op2, op3, initial_state, h, t_mesh, params)",
    doc = """Compute a 3-point correlator of operators 'op1', 'op2' and 'op3' with their time arguments defined on a Cartesian product 't_mesh' x 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_correlator_3t",
    signature = "correlator_3t_container_t(static_operator_t op1, static_operator_t op2, static_operator_t op3,"
                "init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh,"
                "solver_parameters_t params)",
    calling_pattern = "auto result = compute_correlator_3t(op1, op2, op3, initial_state, h, t_mesh, params)",
    doc = """Compute a 3-point correlator of operators 'op1', 'op2' and 'op3' with their time arguments defined on a Cartesian product 't_mesh' x 't_mesh' x 't_mesh'"""
)

#
# compute_g_l()
#

module.add_function(
    name = "compute_g_l",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh, solver_parameters_t params)",
    calling_pattern = "auto result = compute_g_l(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the lesser GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_g_l",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh, solver_parameters_t params)",
    calling_pattern = "auto result = compute_g_l(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the lesser GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

#
# compute_g_g()
#

module.add_function(
    name = "compute_g_g",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh, solver_parameters_t params)",
    calling_pattern = "auto result = compute_g_g(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the greater GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_g_g",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh, solver_parameters_t params)",
    calling_pattern = "auto result = compute_g_g(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the greater GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

#
# compute_g_chi()
#

module.add_function(
    name = "compute_chi",
    signature = "gf_2t_t(chi_indices_t chi_indices, init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh, solver_parameters_t params)",
    calling_pattern = "auto result = compute_chi(chi_indices, initial_state, h, t_mesh, params)",
    doc = """Compute the susceptibility defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_chi",
    signature = "gf_2t_t(chi_indices_t chi_indices, init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh, solver_parameters_t params)",
    calling_pattern = "auto result = compute_chi(chi_indices, initial_state, h, t_mesh, params)",
    doc = """Compute the susceptibility defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

#
# realevol::solver_parameters_t
#

conv = converter_(
    c_type = "realevol::solver_parameters_t ",
    doc = r"""Miscellaneous parameters passed to compute_*() functions"""
)

conv.add_member(c_name = "verbosity",
                c_type = "int",
                initializer = "((mpi::communicator().rank() == 0) ? 3 : 0)",
                doc = r"""Verbosity level""")

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

#
# Comparison tests
#

module.add_function(name = "assert_gfs_are_close",
                    signature = "void(gf_2t_view x, gf_2t_view y, double precision=1.e-6)",
                    doc = """Compare two real time GFs""")

module.add_function(name = "assert_block_gfs_are_close",
                    signature = "void(block_gf_2t_view x, block_gf_2t_view y, double precision=1.e-6)",
                    doc = """Compare two real time block GFs""")

module.generate_code()
