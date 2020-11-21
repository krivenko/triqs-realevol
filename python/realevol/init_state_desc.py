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

# The module
module = module_(full_name = "init_state", app_name = "realevol",
                 doc = "Functions to produce initial states for real-time evolution")

module.add_imports('realevol.operators_texpr', 'realevol.operators_tinterp')

module.add_include("<realevol/init_state.hpp>")

module.add_include("<triqs/cpp2py_converters.hpp>")
module.add_include("<cpp2py/converters/set.hpp>")
module.add_include("<cpp2py/converters/map.hpp>")
module.add_include("<cpp2py/converters/vector.hpp>")
module.add_include("<cpp2py/converters/variant.hpp>")

module.add_preamble("""
using namespace realevol;
using namespace realevol::hilbert_space;
""")

# The class solver
c = class_(
        py_type = "InitState",  # name of the python class
        c_type = "init_state",  # name of the C++ class
        c_type_absolute = "realevol::init_state",
        is_printable = True,
        hdf5 = True,
        doc = "Initial state, including information about the Hilbert space structure"
)

module.add_class(c)

make_pure_init_state_doc = """\
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| Argument Name   | Type                       | Default   | Documentation                                                            |
+=================+============================+===========+==========================================================================+
| generator       | %-26s | --        | Generator operator                                                       |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| fermion_indices | set((str, int))            | --        | Set of (block_index,inner_index) pairs for fermionic degrees of freedom  |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| boson_indices   | set((str, int))            | --        | Set of (block_index,inner_index) pairs for bosonic degrees of freedom    |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| bits_per_boson  | dict((str, int) : int)     | {}        | Number of bits to represent excited states of bosonic degrees of freedom |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
"""

make_equilibrium_init_state_doc = """\
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| Argument Name   | Type                       | Default   | Documentation                                                            |
+=================+============================+===========+==========================================================================+
| h               | %-26s | --        | Hamiltonian in equilibrium                                               |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| fermion_indices | set((str, int))            | --        | Set of (block_index,inner_index) pairs for fermionic degrees of freedom  |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| boson_indices   | set((str, int))            | --        | Set of (block_index,inner_index) pairs for bosonic degrees of freedom    |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| temperature     | float                      | --        | Temperature (zero is a valid value)                                      |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| params          | dict()                     | See below | Equilibrium ED solver parameters                                         |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+
| bits_per_boson  | dict((str, int) : int)     | {}        | Number of bits to represent excited states of bosonic degrees of freedom |
+-----------------+----------------------------+-----------+--------------------------------------------------------------------------+

Equilibrium ED solver parameters
--------------------------------

+------------------------+-----------------+-----------------+------------------------------------------------------------------------------+
| Parameter Name         | Type            | Default         | Documentation                                                                |
+========================+=================+=================+==============================================================================+
| verbosity              | int             | 0               | Verbosity level for the equilibrium solver                                   |
+------------------------+-----------------+-----------------+------------------------------------------------------------------------------+
| min_rel_weight         | float           | Machine epsilon | Discard states with statistical weight below this threshold                  |
+------------------------+-----------------+-----------------+------------------------------------------------------------------------------+
| arpack_min_matrix_size | int             | 101             | Call ARPACK to diagonalize matrices of this size or bigger (cannot be < 4)   |
+------------------------+-----------------+-----------------+------------------------------------------------------------------------------+
| arpack_tolerance       | float           | 0               | Eigenvalue convergence tolerance for ARPACK                                  |
+------------------------+-----------------+-----------------+------------------------------------------------------------------------------+
| arpack_ncv             | dict(int : int) | {}              | ARPACK parameter NCV (number of Lanczos vectors) for each invariant subspace |
+------------------------+-----------------+-----------------+------------------------------------------------------------------------------+
"""

operator_types = [("time_expr_operator_t", "operators_texpr.Operator"),
                  ("time_interp_operator_t", "operators_tinterp.Operator")]

for cpp_t, py_t in operator_types:
    module.add_function(name = "make_pure_init_state",
                        signature = "init_state(%s generator, std::set<indices_t> fermion_indices,"
                                    "std::set<indices_t> boson_indices,"
                                    "std::map<operators::indices_t, int> bits_per_boson)" % cpp_t,
                        calling_pattern = "auto result = make_pure_init_state(generator, fundamental_operator_set(fermion_indices,boson_indices), bits_per_boson)",
                        doc = make_pure_init_state_doc % py_t)

    module.add_function(name = "make_equilibrium_init_state",
                        signature = "init_state(%s h, std::set<indices_t> fermion_indices,"
                                    "std::set<indices_t> boson_indices,"
                                    "double temperature, realevol::eq_solver_parameters_t params,"
                                    "std::map<operators::indices_t, int> bits_per_boson = {})" % cpp_t,
                        calling_pattern = "auto result = make_equilibrium_init_state(h, fundamental_operator_set(fermion_indices,boson_indices),"
                                          "temperature, params, bits_per_boson)",
                        doc = make_equilibrium_init_state_doc % py_t)

conv = converter_(
    c_type = "realevol::eq_solver_parameters_t",
    doc = r"""Parameters of make_equilibrium_init_state()""",
)

conv.add_member(c_name = "verbosity",
                c_type = "int",
                initializer = "0",
                doc = r"""Verbosity level for the equilibrium solver""")

conv.add_member(c_name = "min_rel_weight",
                c_type = "double",
                initializer = "std::numeric_limits<double>::epsilon()",
                doc = r"""Discard states with relative statistical weight below this threshold""")

conv.add_member(c_name = "arpack_min_matrix_size",
                c_type = "int",
                initializer = "101",
                doc = r"""Call ARPACK to diagonalize matrices of this size or bigger (must be >=4)""")

conv.add_member(c_name = "arpack_tolerance",
                c_type = "double",
                initializer = "0",
                doc = r"""Eigenvalue convergence tolerance for ARPACK""")

conv.add_member(c_name = "arpack_ncv",
                c_type = "std::map<long, int>",
                initializer = "{}",
                doc = r"""ARPACK parameter NCV (number of Lanczos vectors) for each invariant subspace""")

module.add_converter(conv)

module.generate_code()
