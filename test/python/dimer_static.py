# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2023, I. Krivenko, M. Danilov, P. Kubiczek
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
import itertools as it
import numpy as np
from scipy.linalg import eigh, expm
from numpy.testing import assert_array_almost_equal

import triqs.utility.mpi as mpi
from triqs.gf import MeshReTime, MeshProduct, Gf, BlockGf

from realevol.operators_tinterp import *
from realevol.init_state import make_equilibrium_init_state
from realevol.realevol import *

#
# Hubbard dimer
#

# Make annihilation matrix operator
def make_c_op(n_fermions, index):
    dim = 2 ** n_fermions
    states = list(it.product((0, 1), repeat = n_fermions))
    c = np.zeros((dim, dim))
    for m, n in it.product(range(dim), range(dim)):
        if states[m] == tuple((states[n][i]-1 if i == index else states[n][i])
                              for i in range(n_fermions)):
            c[m, n] = (-1) ** sum(states[n][:index])
    return c

class dimer_static_ref:

    def make_H(self, beta, U, mu, h, hop):
        # Hamiltonian of decoupled atoms
        H = -mu*(self.n['dn', 0] + self.n['dn', 1])
        H += -mu*(self.n['up', 0] + self.n['up', 1])
        H += U * self.n['dn', 0] @ self.n['up', 0]
        H += U * self.n['dn', 1] @ self.n['up', 1]
        H += 0.5*h*(self.n['up', 0] - self.n['dn', 0])
        H += 0.5*h*(self.n['up', 1] - self.n['dn', 1])
        # Hopping terms
        for s in ('dn', 'up'):
            H += hop * (self.c_dag[s, 0] @ self.c[s, 1] + self.c_dag[s, 1] @ self.c[s, 0])
        return H

    def __init__(self, hbar, beta, U, mu, h, hop0, hop):
        self.hbar = hbar

        # Annihilation operators
        self.c = {('dn', 0) : make_c_op(4, 0),
                  ('dn', 1) : make_c_op(4, 1),
                  ('up', 0) : make_c_op(4, 2),
                  ('up', 1) : make_c_op(4, 3)}
        # Creation and number of particle operators
        self.c_dag = {}
        self.n = {}
        for op in self.c.keys():
            self.c_dag[op] = np.transpose(self.c[op])
            self.n[op] = self.c_dag[op] @ self.c[op]

        H0 = self.make_H(beta, U, mu, h, hop0)
        H = self.make_H(beta, U, mu, h, hop)

        # Calculate \rho_0
        self.rho0 = expm(-beta * H0)
        self.rho0 = self.rho0 / np.trace(self.rho0)

        # Diagonalize H
        self.E, self.V = eigh(H)
        self.Vdag = self.V.conj().T

        # Transform \rho_0 to the eigenbasis of H
        self.rho0 = self.Vdag @ self.rho0 @ self.V

        # Transform c/c^\dagger/n to the eigenbasis of H
        for i in self.c:
            self.c[i] = self.Vdag @ self.c[i] @ self.V
        for i in self.c_dag:
            self.c_dag[i] = self.Vdag @ self.c_dag[i] @ self.V
        for i in self.n:
            self.n[i] = self.Vdag @ self.n[i] @ self.V

    def U(self, t1, t2):
        return np.diag(np.exp(-1j * (t1-t2) * self.E / self.hbar))

    def compute_expvalue(self, t_mesh):
        op = self.c_dag['up', 0] @ self.c['up', 1] + \
             self.c_dag['up', 1] @ self.c['up', 0]
        u = self.U
        res = Gf(mesh = t_mesh, target_shape = [])
        for t in res.mesh:
            res[t] = np.trace(self.rho0 @ u(0, t) @ op @ u(t, 0))
        return res

    def compute_correlator_2t(self, t_mesh):
        op1 = self.c_dag['dn', 0] @ self.c['up', 0] + \
              self.c_dag['dn', 1] @ self.c['up', 1]
        op2 = self.c_dag['up', 0] @ self.c['dn', 0] + \
              self.c_dag['up', 1] @ self.c['dn', 1]
        u = self.U
        corr = Gf(mesh = MeshProduct(t_mesh, t_mesh), target_shape = [])
        for t1, t2 in corr.mesh:
            val = np.trace(self.rho0 @ u(0, t1) @ op1 @ u(t1, t2) @ op2 @ u(t2, 0))
            corr[t1, t2] = val
        return corr

    def compute_correlator_3t(self, t_mesh):
        op1 = self.n['dn', 0] + self.n['up', 0]
        op2 = self.c_dag['dn', 0] @ self.c['up', 0] + \
              self.c_dag['dn', 1] @ self.c['up', 1]
        op3 = self.c_dag['up', 0] @ self.c['dn', 0] + \
              self.c_dag['up', 1] @ self.c['dn', 1]
        u = self.U
        corr = Gf(mesh = MeshProduct(t_mesh, t_mesh, t_mesh), target_shape = [])
        for t1, t2, t3 in corr.mesh:
            val = np.trace(self.rho0 @ u(0, t1) @ op1 @ u(t1, t2) @ op2 @ u(t2, t3) @ op3 @ u(t3, 0))
            corr[t1, t2, t3] = val
        return corr

    def compute_g_g(self, t_mesh):
        u = self.U
        blocks = []
        for s in ('dn', 'up'):
            blocks.append(Gf(mesh = MeshProduct(t_mesh, t_mesh), indices = [0,1]))
            for i, j in it.product((0,1), (0,1)):
                op1 = self.c[s, i]
                op2 = self.c_dag[s, j]
                for t1, t2 in blocks[-1].mesh:
                    val = -1j * np.trace(self.rho0 @ u(0, t1) @ op1 @ u(t1, t2) @ op2 @ u(t2, 0))
                    blocks[-1][t1, t2][i, j] = val / self.hbar
        return BlockGf(name_list = ['dn', 'up'], block_list = blocks)

    def compute_g_l(self, t_mesh):
        u = self.U
        blocks = []
        for s in ('dn', 'up'):
            blocks.append(Gf(mesh = MeshProduct(t_mesh, t_mesh), indices = [0,1]))
            for i, j in it.product((0,1), (0,1)):
                op1 = self.c_dag[s, i]
                op2 = self.c[s, j]
                for t1, t2 in blocks[-1].mesh:
                    val = 1j * np.trace(self.rho0 @ u(0, t2) @ op1 @ u(t2, t1) @ op2 @ u(t1, 0))
                    blocks[-1][t1, t2][i, j] = val / self.hbar
        return BlockGf(name_list = ['dn', 'up'], block_list = blocks)

    def compute_chi(self, t_mesh):
        u = self.U
        chi_indices = list(it.product(['dn', 'up'], [0, 1]))
        chi = Gf(mesh = MeshProduct(t_mesh, t_mesh), indices = list(range(4)))
        for (n1, i1), (n2, i2) in it.product(enumerate(chi_indices), enumerate(chi_indices)):
            op1 = self.n[i1]
            op2 = self.n[i2]
            for t1, t2 in chi.mesh:
                val = -1j * np.trace(self.rho0 @ u(0, t1) @ op1 @ u(t1, t2) @ op2 @ u(t2, 0))
                chi[t1, t2][n1, n2] = val / self.hbar
        return chi

class test_dimer_static(unittest.TestCase):

    hbar = 3.5  # Planck's constant
    beta = 2.0  # Inverse temperature in the initial state
    U = 3.0     # Hubbard interaction
    mu = 1.7    # Chemical potential
    h = 0.2     # Magnetic field
    hop = 0.15  # Hopping amplitude

    @classmethod
    def setUpClass(cls):
        # Fundamental operator set
        cls.fops = set(it.product(('dn', 'up'), (0, 1)))
        cls.gf_struct = [('dn', 2), ('up', 2)]
        cls.chi_indices = list(it.product(['dn', 'up'], [0, 1]))

        # Hamiltonian of decoupled atoms
        cls.H1 = -cls.mu * sum(n(*i) for i in cls.fops)
        for o in (0, 1):
            cls.H1 += cls.U * n('dn', o) * n('up', o)
            cls.H1 += 0.5 * cls.h * (n('up', o) - n('dn', o))

        # H2 = H1 + hopping terms
        cls.H2 = -cls.mu * sum(n(*i) for i in cls.fops)
        for o in (0, 1):
            cls.H2 += cls.U * n('dn', o) * n('up', o)
            cls.H2 += 0.5 * cls.h * (n('up', o) - n('dn', o))
        for s in ('dn', 'up'):
            cls.H2 += cls.hop * (c_dag(s, 0) * c(s, 1) + c_dag(s, 1) * c(s, 0))

        # Thermal initial states generated by H1 and H2
        cls.init_state1 = make_equilibrium_init_state(cls.H1,
                                                fermion_indices = cls.fops,
                                                boson_indices = set(),
                                                temperature = 1 / cls.beta,
                                                params = {'verbosity' : 0})
        cls.init_state2 = make_equilibrium_init_state(cls.H2,
                                                fermion_indices = cls.fops,
                                                boson_indices = set(),
                                                temperature = 1 / cls.beta,
                                                params = {'verbosity' : 0})

        # Real time mesh
        cls.t_mesh = MeshReTime(0, 5.0, 51)

        # Precompute bits needed to produce reference results

        # Hopping switching 0 -> hop
        cls.ref12 = dimer_static_ref(cls.hbar, cls.beta, cls.U, cls.mu, cls.h, 0, cls.hop)
        # Hopping switching hop -> 0
        cls.ref21 = dimer_static_ref(cls.hbar, cls.beta, cls.U, cls.mu, cls.h, cls.hop, 0)

    def test_expvalue(self):
        op = c_dag('up', 0) * c('up', 1) + c_dag('up', 1) * c('up', 0)

        params = {'verbosity' : 0, 'hbar' : self.hbar}

        avg_hop12 = compute_expectval(op, self.init_state1, self.H2, self.t_mesh, params)
        ref12 = self.ref12.compute_expvalue(self.t_mesh)
        assert_array_almost_equal(ref12.data, avg_hop12.data)

        avg_hop21 = compute_expectval(op, self.init_state2, self.H1, self.t_mesh, params)
        ref21 = self.ref21.compute_expvalue(self.t_mesh)
        assert_array_almost_equal(ref21.data, avg_hop21.data)

        # Batch mode
        avg_batch = compute_expectval([2.0 * op, 3.0 * op],
                                      self.init_state1,
                                      self.H2,
                                      self.t_mesh,
                                      params)
        self.assertEqual(len(avg_batch), 2)
        assert_array_almost_equal(2.0 * ref12.data, avg_batch[0].data)
        assert_array_almost_equal(3.0 * ref12.data, avg_batch[1].data)

    def test_compute_correlator_2t(self):
        op1 = c_dag('dn', 0) * c('up', 0) + c_dag('dn', 1) * c('up', 1)
        op2 = c_dag('up', 0) * c('dn', 0) + c_dag('up', 1) * c('dn', 1)

        params = {'verbosity' : 0, 'hbar' : self.hbar}

        corr12 = compute_correlator_2t(op1, op2, self.init_state1, self.H2, self.t_mesh, params)
        ref12 = self.ref12.compute_correlator_2t(self.t_mesh)
        self.assertEqual(ref12.mesh, corr12.mesh)
        self.assertEqual(ref12.target_shape, corr12.target_shape)
        assert_array_almost_equal(ref12.data, corr12.data)

        corr21 = compute_correlator_2t(op1, op2, self.init_state2, self.H1, self.t_mesh, params)
        ref21 = self.ref21.compute_correlator_2t(self.t_mesh)
        self.assertEqual(ref21.mesh, corr21.mesh)
        self.assertEqual(ref21.target_shape, corr21.target_shape)
        assert_array_almost_equal(ref21.data, corr21.data)

        # Batch mode
        corr_batch = compute_correlator_2t([(op1, op2), (2.0 * op1, 3.0 * op2)],
                                           self.init_state1,
                                           self.H2,
                                           self.t_mesh,
                                           params)
        self.assertEqual(len(corr_batch), 2)
        assert_array_almost_equal(ref12.data, corr_batch[0].data)
        assert_array_almost_equal(6.0 * ref12.data, corr_batch[1].data)

    def test_compute_correlator_3t(self):
        op1 = n('dn', 0) + n('up', 0)
        op2 = c_dag('dn', 0) * c('up', 0) + c_dag('dn', 1) * c('up', 1)
        op3 = c_dag('up', 0) * c('dn', 0) + c_dag('up', 1) * c('dn', 1)

        params = {'verbosity' : 0, 'hbar' : self.hbar}

        corr12 = compute_correlator_3t(op1, op2, op3, self.init_state1, self.H2, self.t_mesh, params)
        ref12 = self.ref12.compute_correlator_3t(self.t_mesh)
        self.assertEqual(ref12.mesh, corr12.mesh)
        self.assertEqual(ref12.target_shape, corr12.target_shape)
        assert_array_almost_equal(ref12.data, corr12.data)

        corr21 = compute_correlator_3t(op1, op2, op3, self.init_state2, self.H1, self.t_mesh, params)
        ref21 = self.ref21.compute_correlator_3t(self.t_mesh)
        self.assertEqual(ref21.mesh, corr21.mesh)
        self.assertEqual(ref21.target_shape, corr21.target_shape)
        assert_array_almost_equal(ref21.data, corr21.data)

        # Batch mode
        corr_batch = compute_correlator_3t([(op1, op2, op3), (2.0 * op1, 3.0 * op2, 4.0 * op3)],
                                           self.init_state1,
                                           self.H2,
                                           self.t_mesh,
                                           params)
        self.assertEqual(len(corr_batch), 2)
        assert_array_almost_equal(ref12.data, corr_batch[0].data)
        assert_array_almost_equal(24.0 * ref12.data, corr_batch[1].data)

    def test_g_g(self):
        params = {'verbosity' : 0, 'hbar' : self.hbar}

        g_g12 = compute_g_g(self.gf_struct, self.init_state1, self.H2, self.t_mesh, params)
        ref12 = self.ref12.compute_g_g(self.t_mesh)
        assert_block_gfs_are_close(ref12, g_g12)

        g_g21 = compute_g_g(self.gf_struct, self.init_state2, self.H1, self.t_mesh, params)
        ref21 = self.ref21.compute_g_g(self.t_mesh)
        assert_block_gfs_are_close(ref21, g_g21)

    def test_g_l(self):
        params = {'verbosity' : 0, 'hbar' : self.hbar}

        g_l12 = compute_g_l(self.gf_struct, self.init_state1, self.H2, self.t_mesh, params)
        ref12 = self.ref12.compute_g_l(self.t_mesh)
        assert_block_gfs_are_close(ref12, g_l12)

        g_l21 = compute_g_l(self.gf_struct, self.init_state2, self.H1, self.t_mesh, params)
        ref21 = self.ref21.compute_g_l(self.t_mesh)
        assert_block_gfs_are_close(ref21, g_l21)

    def test_chi(self):
        params = {'verbosity' : 0, 'hbar' : self.hbar}

        chi12 = compute_chi(self.chi_indices, self.init_state1, self.H2, self.t_mesh, params)
        ref12 = self.ref12.compute_chi(self.t_mesh)
        assert_gfs_are_close(ref12, chi12)

        chi21 = compute_chi(self.chi_indices, self.init_state2, self.H1, self.t_mesh, params)
        ref21 = self.ref21.compute_chi(self.t_mesh)
        assert_gfs_are_close(ref21, chi21)

if __name__ == '__main__':
    unittest.main()
