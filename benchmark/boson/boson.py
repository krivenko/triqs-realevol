from realevol.texpr import TExpr
from realevol.operators_texpr import c, c_dag, n, a, a_dag
from realevol.init_state import *
from realevol import Solver

gf_struct = {'up':[0], 'dn':[0]}
chi_indices = [('dn',0),('up',0)]
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
S.initial_state = init_state

S.compute_2t_obs(h = h, params = {})
