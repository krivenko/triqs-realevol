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
from triqs.gf import BlockGf, MeshReTime
from realevol.tinterp import TInterp as ti
from realevol.operators_tinterp import *
from realevol.init_state import *
from realevol.realevol import *
from itertools import product
import numpy as np

# Cluster benchmark
# 4 site cluster, nn-hoppings are switched on (quench and slow switching)

spin_names = ('up','dn')

t_max = 1.0
n_t = 11
m_interp = MeshReTime(0, t_max, n_t)

# Model parameters
U = 2.0
mu = U*0.1
h = 0.15
t = 0.3
tp = ti(m_interp, np.array([-0.1*(1-np.exp(-10*x)) for x in m_interp]))
beta = 40.0

gf_struct = {'dn' : range(3)}
chi_indices = [(sn,site) for sn, site in product(spin_names,range(4))]

fops = set((sn,site) for sn, site in product(spin_names,range(4)))
print("Fundamental operator set:", fops)

eq_params = {}
eq_params['verbosity'] = 2
eq_params['arpack_min_matrix_size'] = 11

# Initial Hamiltonian
h0 = mu*sum(n(*i) for i in fops)
h0 += sum((h * n(*i) if i[0] == "up" else -h * n(*i)) for i in fops)
h0 += U*sum(n('up',site)*n('dn',site) for site in range(4))
h0 += t*sum(c_dag(sn,site)*c(sn,(site+1)%4) + c_dag(sn,site)*c(sn,(site-1)%4)
            for sn, site in product(spin_names,range(4)))

print("h0 =", h0)

init_state = make_equilibrium_init_state(h0, fermion_indices = fops, boson_indices = set(),
                                         temperature = 1/beta, params = eq_params)

print("Initial state:")
print(init_state)

# Time-dependent Hamiltonian
h = h0 + tp*sum(c_dag(sn,site1)*c(sn,site2)
                for sn, (site1, site2) in product(spin_names,((0,2),(2,0),(1,3),(3,1))))

print("h(t) =", h)

# Solver object
S = Solver(gf_struct, chi_indices, t_max = t_max, n_t = n_t)

# Set initial state
S.set_initial_state(init_state)

gf_params = {}
gf_params['verbosity'] = 2
gf_params['lanczos_min_matrix_size'] = 10000

S.compute_2t_obs(h = h, params = gf_params)

if mpi.is_master_node():
    with HDFArchive('cluster.h5', 'w') as ar:
        ar['init_state'] = init_state
        ar['g_l'] = S.g_l
        ar['g_g'] = S.g_g
        ar['chi'] = S.chi
