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

import unittest

from h5 import HDFArchive
from triqs.gf import BlockGf, MeshReTime
from realevol.tinterp import TInterp as ti
from realevol.operators_tinterp import *
from realevol.init_state import *
from realevol.realevol import *
import triqs.utility.mpi as mpi
import numpy as np
from itertools import product

class test_star_texpr(unittest.TestCase):
    """Star test: correlated Hubbard atom + n bath sites"""

    def test_Solver(self):
        spin_names = ('up','dn')

        # Fine mesh to construct interpolation objects
        m_interp = MeshReTime(0, 5.0, 1001)

        # Model parameters
        U = 3.0
        mu = U*0.5
        n_bath = 2
        eps = np.linspace(-0.2,0.2,n_bath)
        t = [0.3]*n_bath
        dt = ti(m_interp, np.array([0.1*(1-np.exp(-5*x)) for x in m_interp]))

        fops = set(product(spin_names,range(n_bath+1)))
        gf_struct = {'dn':[0], 'up':[0]}
        chi_indices = [('dn',0),('up',0)]

        t_max = 5.0
        n_t = 6

        ## Initial Hamiltonian
        h0 = -mu*(n('up',0) + n('dn',0)) + U*n('up',0)*n('dn',0)
        h0 += sum(eps[i-1]*(n('up',i)+n('dn',i)) for i in range(1,n_bath+1))
        h0 += sum(-t[i-1]*(c_dag(sn,0)*c(sn,i) + c_dag(sn,i)*c(sn,0))
                  for sn, i in product(spin_names,range(1,n_bath+1)))

        init_state = make_equilibrium_init_state(h0,
                                                fermion_indices = fops,
                                                boson_indices = set(),
                                                temperature = 0,
                                                params = {'verbosity':2})

        # Hamiltonian after quench
        h = h0 + sum(dt*(c_dag(sn,0)*c(sn,1) + c_dag(sn,1)*c(sn,0)) for sn in spin_names)

        # Solver object
        S = Solver(gf_struct, chi_indices, t_max = t_max, n_t = n_t)

        # Set initial state
        S.set_initial_state(init_state)

        gf_params = {}
        gf_params['verbosity'] = 2
        gf_params['lanczos_min_matrix_size'] = 10000

        S.compute_2t_obs(h = h, params = gf_params)

        if mpi.is_master_node():
            with HDFArchive('star.ref.h5','r') as ar:
                #ar['init_state'] = init_state
                #ar['g_l'] = S.g_l
                #ar['g_g'] = S.g_g
                #ar['chi'] = S.chi
                assert_block_gfs_are_close(ar['g_g'], S.g_g)
                assert_block_gfs_are_close(ar['g_l'], S.g_l)
                assert_gfs_are_close(ar['chi'], S.chi)

if __name__ == '__main__':
    unittest.main()
