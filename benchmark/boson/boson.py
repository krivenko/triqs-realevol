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

from h5 import HDFArchive
import triqs.utility.mpi as mpi
from realevol.texpr import TExpr
from realevol.operators_texpr import c, c_dag, n, a, a_dag
from realevol.init_state import *
from realevol.realevol import Solver

gf_struct = [('dn', [0]), ('up', [0])]
chi_indices = [('dn', 0),('up',0)]
fops_fermion = set([("up",0),("dn",0)])
fops_boson = set([("B",0)])
bits_per_boson = {('B',0) : 3}
t_max = 10.0
n_t = 1000

S = Solver(gf_struct, chi_indices, t_max = t_max, n_t = n_t)

U = 3.0
mu = U/2
Omega = 0.5
Lambda = TExpr("0.2*sgn(t+1)")

h0 = -mu*(n("up",0) + n("dn",0)) + U*n("up",0)*n("dn",0)
h0 += Omega*a_dag("B",0)*a("B",0)

h = -mu*(n("up",0) + n("dn",0)) + U*n("up",0)*n("dn",0)
h += Omega*a_dag("B",0)*a("B",0)
h += Lambda*(n("up",0) + n("dn",0))*(a_dag("B",0) + a("B",0))

init_state = make_equilibrium_init_state(h0,
                                         fermion_indices = fops_fermion,
                                         boson_indices = fops_boson,
                                         temperature = 0,
                                         bits_per_boson = bits_per_boson,
                                         params = {})

# Set initial state
S.set_initial_state(init_state)

S.compute_2t_obs(h = h, params = {})

if mpi.is_master_node():
    with HDFArchive('boson.h5', 'w') as ar:
        ar['init_state'] = init_state
        ar['g_l'] = S.g_l
        ar['g_g'] = S.g_g
        ar['chi'] = S.chi
