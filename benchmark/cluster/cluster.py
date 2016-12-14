from realevol.texpr import TExpr as te
from realevol.operators import *
from realevol.init_state import *
from itertools import product

# Cluster benchmark
# 4 site cluster, nn-hoppings are switched on (quench and slow switching)

spin_names = ('up','dn')

# Model parameters
U = 2.0
mu = U*0.1
t = 0.3
tp = -0.1
beta = 20.0

fops = set((sn,site) for sn, site in product(spin_names,range(4)))
print "Fundamental operator set:", fops

eq_params = {}
eq_params['verbosity'] = 2
eq_params['arpack_min_matrix_size'] = 11

# Initial Hamiltonian
h0 = mu*sum(n(*i) for i in fops)
h0 += U*sum(n('up',site)*n('dn',site) for site in range(4))
h0 += t*sum(c_dag(sn,site)*c(sn,(site+1)%4) + c_dag(sn,site)*c(sn,(site-1)%4)
            for sn, site in product(spin_names,range(4)))

print "h0 =", h0

init_state = make_equilibrium_init_state(h0, fermion_indices = fops, boson_indices = set(),
                                         temperature = 1/beta, params = eq_params)

print "Initial state:"
print init_state
