from realevol.texpr import TExpr
from realevol.operators import c, c_dag, n, a, a_dag
from realevol import Solver

gf_struct = {'up':[0], 'dn':[0]}
time_window = (0.0,10.0)
n_t = 1000

S = Solver(gf_struct, time_window = time_window, n_t = n_t)

U = 3.0
mu = U/2
Omega = 0.5
Lambda = TExpr("0.2*sgn(t+1)")

h = -mu*(n("up",0) + n("dn",0)) + U*n("up",0)*n("dn",0)
h += Omega*a_dag("B",0)*a("B",0)
h += Lambda*(n("up",0) + n("dn",0))*(a_dag("B",0) + a("B",0))

h0 = -mu*(n("up",0) + n("dn",0)) + U*n("up",0)*n("dn",0)
h0 += Omega*a_dag("B",0)*a("B",0)

params = {}
params['h'] = h
params['h0'] = h0
params['bits_per_boson'] = {('B',0) : 3}

print h
S.solve(**params)
