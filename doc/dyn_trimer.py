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

import numpy as np
from scipy.integrate import complex_ode
from itertools import product

from h5 import *

M = 6
fock_dim = 2 ** M

states = [s for s in product(*[range(0,2)]*M)]
# Sort by the total number of particles
states.sort(key = lambda s: sum(s))
states = list(map(lambda s: np.array(s,dtype=int), states))

ct = [np.matrix(np.zeros((fock_dim,fock_dim)),dtype=complex) for o in range(M)]
cdt = [np.matrix(np.zeros((fock_dim,fock_dim)),dtype=complex) for o in range(M)]

for o in range(M):
    mask = np.zeros(M)
    mask[o] = 1
    for m,n in product(range(fock_dim),range(fock_dim)):
        if (states[n] - states[m] == mask).all(): ct[o][m,n] = (-1.0) ** sum(states[n][0:o])
    cdt[o] = np.transpose(ct[o])

canonical_check = np.zeros((M,M),dtype=bool)
for o1, o2 in product(range(M),range(M)):
    canonical_check[o1,o2] = (cdt[o1]*ct[o2] + ct[o2]*cdt[o1] == np.eye(fock_dim)*(1 if o1==o2 else 0)).all()

print("Canonical anticommutators:")
print(canonical_check)

c1up, c2up, c3up, c1dn, c2dn, c3dn = ct
cd1up, cd2up, cd3up, cd1dn, cd2dn, cd3dn = cdt

n1up, n2up, n3up = cd1up*c1up, cd2up*c2up, cd3up*c3up
n1dn, n2dn, n3dn = cd1dn*c1dn, cd2dn*c2dn, cd3dn*c3dn

hbar = 1.0
U = 1.0
mu = U/2
V = lambda t: 0.3*(1.0 - np.exp(-4.0*t))

tmin, tmax = 0.0, 50.0
N_t = 1001
dt = (tmax-tmin)/(N_t-1)

def H(t):
    tmp = -mu*(n1up + n2up + n3up + n1dn + n2dn + n3dn)
    tmp += U*(n1up*n1dn + n2up*n2dn + n3up*n3dn)
    tmp += -V(t)*(cd1up*c2up + cd2up*c1up + cd1up*c3up + cd3up*c1up + cd2up*c3up + cd3up*c2up)
    tmp += -V(t)*(cd1dn*c2dn + cd2dn*c1dn + cd1dn*c3dn + cd3dn*c1dn + cd2dn*c3dn + cd3dn*c2dn)
    return tmp

SO = complex_ode(lambda t,y: np.dot(H(t),y)/(1j*hbar))
SO

vac = np.zeros(fock_dim,dtype=complex)
vac[0] = 1.0

def average(op,psi):
    res = np.zeros(len(psi))
    for n, psi_at_t in enumerate(psi):
        res[n] = np.einsum("i,ij,j",np.conjugate(psi_at_t),op,psi_at_t).real
    return res

arch = HDFArchive('trimer.output.h5','w')

# S_z = 1/2
psi0 = np.dot(cd1up*cd2up*cd3dn,vac)
SO.set_initial_value(psi0,t=tmin)

psi = [np.ravel(psi0)]
while SO.successful() and SO.t <= tmax-dt:
    SO.integrate(SO.t+dt)
    psi.append(SO.y)

avg_names = ("n_1_up","n_2_up","n_3_up","n_1_dn","n_2_dn","n_3_dn","unity")
avgs = list(map(lambda op: average(op,psi), (n1up,n2up,n3up,n1dn,n2dn,n3dn)))
avgs.append(np.ones(N_t))
arch.create_group("Sz_0.5")
Sz_gr = arch['Sz_0.5']

for name, avg in zip(avg_names,avgs):
    Sz_gr.create_group(name)
    Sz_gr[name]['vector'] = avg
    Sz_gr[name].create_group('mesh')
    Sz_gr[name]['mesh']['min'] = tmin
    Sz_gr[name]['mesh']['max'] = tmax
    Sz_gr[name]['mesh']['size'] = N_t

# S_z = 0
psi0 = np.dot(cd1dn*cd1up*cd2up*cd3dn,vac)
SO.set_initial_value(psi0,t=tmin)

psi = [np.ravel(psi0)]
while SO.successful() and SO.t <= tmax-dt:
    SO.integrate(SO.t+dt)
    psi.append(SO.y)

avg_names = ("n_1_up","n_2_up","n_3_up","n_1_dn","n_2_dn","n_3_dn","unity")
avgs = list(map(lambda op: average(op,psi), (n1up,n2up,n3up,n1dn,n2dn,n3dn)))
avgs.append(np.ones(N_t))
arch.create_group("Sz_0")
Sz_gr = arch['Sz_0']

for name, avg in zip(avg_names,avgs):
    Sz_gr.create_group(name)
    Sz_gr[name]['vector'] = avg
    Sz_gr[name].create_group('mesh')
    Sz_gr[name]['mesh']['min'] = tmin
    Sz_gr[name]['mesh']['max'] = tmax
    Sz_gr[name]['mesh']['size'] = N_t
