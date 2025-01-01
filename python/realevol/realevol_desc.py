# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2025, I. Krivenko, M. Danilov, P. Kubiczek
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
module.add_include("<realevol/make_static_op.hpp>")
module.add_include("<realevol/time_expr.hpp>")
module.add_include("<realevol/time_interp.hpp>")
module.add_include("<realevol/utility.hpp>")
module.add_include("<triqs/gfs/gf_tests.hpp>")

module.add_include("<triqs/cpp2py_converters.hpp>")
module.add_include("<cpp2py/converters/std_array.hpp>")
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

template<typename OperatorType>
std::vector<static_operator_t> make_static_ops(std::vector<OperatorType> const& ops) {
  std::vector<static_operator_t> result;
  result.reserve(ops.size());
  for(auto const& op : ops)
    result.emplace_back(make_static_op(op));
  return result;
}

template<typename OperatorType, std::size_t NPoints>
std::vector<std::array<static_operator_t, NPoints>>
make_static_ops(std::vector<std::array<OperatorType, NPoints>> const& ops) {
  std::vector<std::array<static_operator_t, NPoints>> result;
  result.reserve(ops.size());
  for(auto const& op_array : ops) {
    result.emplace_back(
      map_array<static_operator_t>([](auto const& op) {
        return make_static_op(op);
      }, op_array)
    );
  }
  return result;
}
""")

module.add_enum("h_interpolation", ["Rectangle", "Trapezoid", "Simpson"], "realevol",
                "Hamiltonian interpolation between time slices")

operator_types = [("time_expr_operator_t", "operators_texpr.Operator"),
                  ("time_interp_operator_t", "operators_tinterp.Operator")]

#
# compute_expectval()
#

for cpp_t, py_t in operator_types:

    module.add_function(
        name = "compute_expectval",
        signature = f"expectval_container_t({cpp_t} op,"
                    f"init_state initial_state, {cpp_t} h, mesh_t_t t_mesh,"
                    "solver_parameters_t<1> params)",
        calling_pattern = "auto result = compute_expectval(make_static_op(op), initial_state, h, t_mesh, params)",
        doc = """Compute expectation value of operator 'op' as a function of time"""
    )

    # Batch version
    module.add_function(
        name = "compute_expectval",
        signature = f"std::vector<expectval_container_t>(std::vector<{cpp_t}> ops,"
                    f"init_state initial_state, {cpp_t} h, mesh_t_t t_mesh,"
                    "solver_parameters_t<1> params)",
        calling_pattern = "auto result = compute_expectval(make_static_ops(ops), initial_state, h, t_mesh, params)",
        doc = """Batch-compute expectation values of operators 'ops' as functions of time"""
    )

#
# compute_correlator_2t()
#

for cpp_t, py_t in operator_types:

    module.add_function(
        name = "compute_correlator_2t",
        signature = f"correlator_2t_container_t({cpp_t} op1, {cpp_t} op2,"
                    f"init_state initial_state, {cpp_t} h, mesh_t_t t_mesh,"
                    "solver_parameters_t<2> params)",
        calling_pattern = "auto result = compute_correlator_2t(make_static_op(op1), make_static_op(op2), initial_state, h, t_mesh, params)",
        doc = """Compute a 2-point correlator of operators 'op1' and 'op2' with their time arguments defined on a Cartesian product 't_mesh' x 't_mesh'"""
    )

    # Batch version
    module.add_function(
        name = "compute_correlator_2t",
        signature = f"std::vector<correlator_2t_container_t>(std::vector<std::array<{cpp_t}, 2>> ops,"
                    f"init_state initial_state, {cpp_t} h, mesh_t_t t_mesh,"
                    "solver_parameters_t<2> params)",
        calling_pattern = "auto result = compute_correlator_2t(make_static_ops(ops), initial_state, h, t_mesh, params)",
        doc = """Batch-compute 2-point correlators of operators 'ops[i][0]' and 'ops[i][1]' with their time arguments defined """
              """on a Cartesian product 't_mesh' x 't_mesh'"""
    )

#
# compute_correlator_3t()
#

for cpp_t, py_t in operator_types:

    module.add_function(
        name = "compute_correlator_3t",
        signature = f"correlator_3t_container_t({cpp_t} op1, {cpp_t} op2, {cpp_t} op3,"
                    f"init_state initial_state, {cpp_t} h, mesh_t_t t_mesh,"
                    "solver_parameters_t<3> params)",
        calling_pattern = "auto result = compute_correlator_3t(make_static_op(op1), make_static_op(op2), make_static_op(op3), initial_state, h, t_mesh, params)",
        doc = """Compute a 3-point correlator of operators 'op1', 'op2' and 'op3' with their time arguments defined on a Cartesian product 't_mesh' x 't_mesh' x 't_mesh'"""
    )

    # Batch version
    module.add_function(
        name = "compute_correlator_3t",
        signature = f"std::vector<correlator_3t_container_t>(std::vector<std::array<{cpp_t}, 3>> ops,"
                    f"init_state initial_state, {cpp_t} h, mesh_t_t t_mesh,"
                    "solver_parameters_t<3> params)",
        calling_pattern = "auto result = compute_correlator_3t(make_static_ops(ops), initial_state, h, t_mesh, params)",
        doc = """Batch-compute 3-point correlators of operators 'ops[i][0]', 'ops[i][1]' and 'ops[i][2]' """
              """with their time arguments defined on a Cartesian product 't_mesh' x 't_mesh' x 't_mesh'"""
    )

#
# compute_g_l()
#

module.add_function(
    name = "compute_g_l",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh, solver_parameters_t<2> params)",
    calling_pattern = "auto result = compute_g_l(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the lesser GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_g_l",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh, solver_parameters_t<2> params)",
    calling_pattern = "auto result = compute_g_l(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the lesser GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

#
# compute_g_g()
#

module.add_function(
    name = "compute_g_g",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh, solver_parameters_t<2> params)",
    calling_pattern = "auto result = compute_g_g(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the greater GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_g_g",
    signature = "block_gf_2t_t(gf_struct_t gf_struct, init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh, solver_parameters_t<2> params)",
    calling_pattern = "auto result = compute_g_g(gf_struct, initial_state, h, t_mesh, params)",
    doc = """Compute the greater GF defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

#
# compute_g_chi()
#

module.add_function(
    name = "compute_chi",
    signature = "gf_2t_t(chi_indices_t chi_indices, init_state initial_state, time_expr_operator_t h, mesh_t_t t_mesh, solver_parameters_t<2> params)",
    calling_pattern = "auto result = compute_chi(chi_indices, initial_state, h, t_mesh, params)",
    doc = """Compute the susceptibility defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

module.add_function(
    name = "compute_chi",
    signature = "gf_2t_t(chi_indices_t chi_indices, init_state initial_state, time_interp_operator_t h, mesh_t_t t_mesh, solver_parameters_t<2> params)",
    calling_pattern = "auto result = compute_chi(chi_indices, initial_state, h, t_mesh, params)",
    doc = """Compute the susceptibility defined on a Cartesian product 't_mesh' x 't_mesh'"""
)

#
# realevol::solver_parameters_t
#

for npoints in range(1, 4):
    conv = converter_(
        c_type = f"realevol::solver_parameters_t<%d> " % npoints,
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

    conv.add_member(c_name = "t_ranges",
                    c_type = f"std::array<std::pair<double, double>, %i>" % npoints,
                    initializer = f"make_array_repeat<std::pair<double, double>, %i>(std::make_pair(-INFINITY, INFINITY))" % npoints,
                    doc = r"""Time argument range restrictions, one per correlator point""")

    conv.add_member(c_name = "delta_t_max",
                    c_type = f"std::array<double, %i>" % (npoints - 1),
                    initializer = f"make_array_repeat<double, %i>(INFINITY)" % (npoints - 1),
                    doc = r"""Maximal separations of successive time arguments""")

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
# Utility functions
#

module.add_function(
    name = "make_gf_ret_adv",
    signature = "std::pair<block_gf_2t_t,block_gf_2t_t>(block_gf_2t_t g_l, block_gf_2t_t g_g)",
    doc = """Compute retarded and advanced Green's functions out of the lesser and greater components"""
)

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
