from pytriqs.gf import BlockGf
from pytriqs.archive import HDFArchive
from realevol.texpr import TExpr as te
from realevol.operators import *
from realevol.init_state import *
from realevol.realevol import *
import pytriqs.utility.mpi as mpi
import numpy as np
from datetime import datetime
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

t_max = 1.0
n_t = 11

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

print "Computing with LAPACK ..."
timing = datetime.now()
S.compute_2t_obs(h = h, **gf_params)
timing = datetime.now() - timing
print "Elapsed time:", timing

lapack = {'g_l' : S.g_l.copy(),
          'g_g' : S.g_g.copy(),
          'chi' : S.chi.copy(),
          'timing' : str(timing)}

print "Computing with Lanczos algorithm ..."
gf_params['lanczos_min_matrix_size'] = 9
timing = datetime.now()
S.compute_2t_obs(h = h, **gf_params)
timing = datetime.now() - timing
print "Elapsed time:", timing

lanczos = {'g_l' : S.g_l.copy(),
           'g_g' : S.g_g.copy(),
           'chi' : S.chi.copy(),
           'timing' : str(timing)}

if mpi.is_master_node():
    with HDFArchive('lapack_vs_lanczos.n_t_%i.h5' % n_t, 'w') as ar:
        ar['init_state'] = init_state
        ar['lapack'] = lapack
        ar['lanczos'] = lanczos
        ar.create_group('diffs')
        gr = ar['diffs']

        diff = lambda x, y: np.max(np.abs(x.data - y.data))
        gr['g_l_up'] = diff(lapack['g_l']['up'], lanczos['g_l']['up'])
        gr['g_l_dn'] = diff(lapack['g_l']['dn'], lanczos['g_l']['dn'])
        gr['g_g_up'] = diff(lapack['g_g']['up'], lanczos['g_g']['up'])
        gr['g_g_dn'] = diff(lapack['g_g']['dn'], lanczos['g_g']['dn'])
        for i, j in product(range(2), range(2)):
            gr['chi_%i_%i' % (i,j)] = diff(lapack ['chi'][i, j],
                                           lanczos['chi'][i, j])
