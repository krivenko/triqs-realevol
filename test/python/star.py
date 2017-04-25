from pytriqs.gf import BlockGf
from pytriqs.archive import HDFArchive
from realevol.texpr import TExpr as te
from realevol.operators import *
from realevol.init_state import *
from realevol.realevol import *
import pytriqs.utility.mpi as mpi
import numpy as np
from itertools import product

# Star test: correlated Hubbard atom + n bath sites

spin_names = ('up','dn')

# Model parameters
U = 3.0
mu = U*0.5
n_bath = 2
eps = np.linspace(-0.2,0.2,n_bath)
t = [0.3]*n_bath
dt = te("0.1*(1-exp(-5*t))")

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
S.initial_state = init_state

gf_params = {}
gf_params['verbosity'] = 2
gf_params['lanczos_min_matrix_size'] = 10000

S.compute_2t_obs(h = h, **gf_params)

if mpi.is_master_node():
    with HDFArchive('star.ref.h5','r') as ar:
        #ar['init_state'] = init_state
        #ar['g_l'] = S.g_l
        #ar['g_g'] = S.g_g
        #ar['chi'] = S.chi
        # FIXME: Workaround for library issue #342
        for bn in spin_names:
            assert_gfs_are_close(ar['g_l'][bn], S.g_l[bn])
            assert_gfs_are_close(ar['g_g'][bn], S.g_g[bn])
        assert_gfs_are_close(ar['chi'], S.chi)

