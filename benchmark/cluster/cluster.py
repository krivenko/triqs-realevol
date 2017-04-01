from realevol.texpr import TExpr as te
from realevol.operators import *
from realevol.init_state import *
from realevol.realevol import *
from itertools import product

# Cluster benchmark
# 4 site cluster, nn-hoppings are switched on (quench and slow switching)

spin_names = ('up','dn')

# Model parameters
U = 2.0
mu = U*0.1
t = 0.3
#tp = -0.1
tp = -0.1*te("1-exp(-10*t)")
beta = 40.0

gf_struct = {'dn' : range(3)}
chi_indices = [(sn,site) for sn, site in product(spin_names,range(4))]
t_max = 1.0
n_t = 11

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

# Time-dependent Hamiltonian
h = h0 + tp*sum(c_dag(sn,site1)*c(sn,site2)
                for sn, (site1, site2) in product(spin_names,((0,2),(2,0),(1,3),(3,1))))

print "h(t) =", h

# Solver object
S = Solver(gf_struct, chi_indices, t_max = t_max, n_t = n_t)

# Set initial state
S.initial_state = init_state

gf_params = {}
gf_params['verbosity'] = 2
gf_params['lanczos_min_matrix_size'] = 10000

S.compute_2t_obs(h = h, **gf_params)
