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

from realevol.tinterp import TInterp, is_constant, is_zero, conj
from triqs.gf import MeshReTime
from numpy import array, where
from numpy.testing import assert_allclose

class test_tinterp(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        m = MeshReTime(0, 10, 11)
        cls.m = m

        cls.TI = [TInterp(),
                  TInterp(m, array([(t*t).real for t in m])),
                  TInterp(m, array([(t + 1).real for t in m])),
                  TInterp(m, 4.5),
                  TInterp(m, array([(t*t).real for t in m]), 1.0),
                  TInterp(m, array([(t + 1).real for t in m]), array([(t*t).real for t in m])),
                  TInterp(m, 4.5, 2.9),
                  TInterp(m, 1.9+2.8j),
                  TInterp(m, array([t*t + (t + 1)*1j for t in m]))]

        cls.T = array([0, 0.1, 1, 5.5, 10])
        cls.TI_ref = [array([0, 0, 0, 0, 0]),
                      array([0, 0.1, 1, 30.5, 100]),
                      array([1, 1.1, 2, 6.5, 11]),
                      array([4.5, 4.5, 4.5, 4.5, 4.5]),
                      array([1.0j, 0.1+1.0j, 1+1.0j, 30.5+1.0j, 100+1.0j]),
                      array([1, 1.1+0.1j, 2+1j, 6.5+30.5j, 11+100j]),
                      array([4.5+2.9j, 4.5+2.9j, 4.5+2.9j, 4.5+2.9j, 4.5+2.9j]),
                      array([1.9+2.8j, 1.9+2.8j, 1.9+2.8j, 1.9+2.8j, 1.9+2.8j]),
                      array([1j, 0.1+1.1j, 1+2j, 30.5+6.5j, 100+11j])]

    @classmethod
    def assertEvalAllClose(cls, ti, ti_ref, tol = 1e-12):
        assert_allclose([ti(t) for t in cls.T], ti_ref, atol = tol)

    def test_is_constant(self):
        self.assertEqual([is_constant(t) for t in self.TI],
                         [True, False, False, True, False, False, True, True, False])

    def test_is_zero(self):
        self.assertTrue(is_zero(self.TI[0]))
        self.assertFalse(is_zero(self.TI[1]))

    def test_conj(self):
        m = self.m
        self.assertEqual(conj(self.TI[3]), TInterp(m, 4.5))
        self.assertEqual(conj(self.TI[5]), TInterp(m, array([t + 1 - t*t*1j for t in m])))

    def test_eval(self):
        for i in range(len(self.TI)):
            self.assertEvalAllClose(self.TI[i], self.TI_ref[i])

    def test_minus(self):
        self.assertEvalAllClose(-self.TI[2], -self.TI_ref[2])
        self.assertEvalAllClose(-self.TI[5], -self.TI_ref[5])

    def test_addition(self):
        self.assertEvalAllClose(self.TI[1] + self.TI[2], self.TI_ref[1] + self.TI_ref[2])
        self.assertEvalAllClose(self.TI[1] + 0.5, self.TI_ref[1] + 0.5)
        self.assertEvalAllClose(0.5 + self.TI[2], 0.5 + self.TI_ref[2])
        self.assertEvalAllClose(self.TI[1] + 0.5j, self.TI_ref[1] + 0.5j)
        self.assertEvalAllClose(0.5j + self.TI[2], 0.5j + self.TI_ref[2])

    def test_subtraction(self):
        self.assertEvalAllClose(self.TI[1] - self.TI[2], self.TI_ref[1] - self.TI_ref[2])
        self.assertEvalAllClose(self.TI[1] - 0.5, self.TI_ref[1] - 0.5)
        self.assertEvalAllClose(0.5 - self.TI[2], 0.5 - self.TI_ref[2])
        self.assertEvalAllClose(self.TI[1] - 0.5j, self.TI_ref[1] - 0.5j)
        self.assertEvalAllClose(self.TI[1] - 0.5j, self.TI_ref[1] - 0.5j)

    def test_multiplication(self):
        self.assertEvalAllClose(self.TI[1] * self.TI[2], array([0, 0.2, 2, 201, 1100]))
        self.assertEvalAllClose(self.TI[5] * self.TI[6],
                                array([4.5+2.9j, 4.66+3.64j, 6.1+10.3j, -59.2+156.1j, -240.5+481.9j]))
        self.assertEvalAllClose(self.TI[1] * 0.5, self.TI_ref[1] * 0.5)
        self.assertEvalAllClose(0.5 * self.TI[2], 0.5 * self.TI_ref[2])
        self.assertEvalAllClose(self.TI[1] * 0.5j, self.TI_ref[1] * 0.5j)
        self.assertEvalAllClose(0.5j * self.TI[2], 0.5j * self.TI_ref[2])

    def test_division(self):
        self.assertEvalAllClose(self.TI[1] / 0.5, self.TI_ref[1] / 0.5)
        self.assertEvalAllClose(self.TI[1] / 0.5j, self.TI_ref[1] / 0.5j)

if __name__ == '__main__':
    unittest.main()
