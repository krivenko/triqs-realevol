# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2022, I. Krivenko, M. Danilov, P. Kubiczek
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

from realevol.texpr import TExpr, is_constant, is_zero, conj
from math import sin, sqrt, pi
from numpy import array
from numpy.testing import assert_equal

class test_texpr(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.T = [0, 0.25, 16, 64]

        cls.TE = [TExpr("t^2"),
                  TExpr("t + sin(pi/2)"),
                  TExpr("sqrt(9.0) + 1.5"),
                  TExpr("t^2",1.0),
                  TExpr("t + sin(pi/2)","t^3"),
                  TExpr("sqrt(9.0) + 1.5","2.9"),
                  TExpr(1.9+2.8j)]

        cls.TE_ref = [lambda t: t**2,
                      lambda t: t+sin(pi/2),
                      lambda t: sqrt(9.0)+1.5,
                      lambda t: t**2+1j,
                      lambda t: t+sin(pi/2)+((t**3)*1j),
                      lambda t: sqrt(9.0)+1.5+2.9j,
                      lambda t: 1.9+2.8j]

    @classmethod
    def vals(cls, f):
        return array([f(t) for t in cls.T], dtype=complex)

    @classmethod
    def assertEvalAllEqual(cls, te, ref):
        assert_equal(cls.vals(te), ref)

    def test_is_constant(self):
        self.assertEqual([is_constant(t) for t in self.TE],
                         [False, False, True, False, False, True, True])

    def test_eval(self):
        for te, te_ref in zip(self.TE, self.TE_ref):
            self.assertEvalAllEqual(te, self.vals(te_ref))

    def test_minus(self):
        self.assertEvalAllEqual(-self.TE[1], -self.vals(self.TE_ref[1]))
        self.assertEvalAllEqual(-self.TE[4], -self.vals(self.TE_ref[4]))

    def test_addition(self):
        self.assertEvalAllEqual(self.TE[0] + self.TE[1],
                                self.vals(self.TE_ref[0]) + self.vals(self.TE_ref[1]))
        self.assertEvalAllEqual(self.TE[0] + 0.5, self.vals(self.TE_ref[0]) + 0.5)
        self.assertEvalAllEqual(0.5 + self.TE[1], 0.5 + self.vals(self.TE_ref[1]))
        self.assertEvalAllEqual(self.TE[0] + 0.5j, self.vals(self.TE_ref[0]) + 0.5j)
        self.assertEvalAllEqual(0.5j + self.TE[1], 0.5j + self.vals(self.TE_ref[1]))

    def test_subtraction(self):
        self.assertEvalAllEqual(self.TE[0] - self.TE[1],
                                self.vals(self.TE_ref[0]) - self.vals(self.TE_ref[1]))
        self.assertEvalAllEqual(self.TE[0] - 0.5, self.vals(self.TE_ref[0]) - 0.5)
        self.assertEvalAllEqual(0.5 - self.TE[1], 0.5 - self.vals(self.TE_ref[1]))
        self.assertEvalAllEqual(self.TE[0] - 0.5j, self.vals(self.TE_ref[0]) - 0.5j)
        self.assertEvalAllEqual(0.5j - self.TE[1], 0.5j - self.vals(self.TE_ref[1]))

    def test_multiplication(self):
        self.assertEvalAllEqual(self.TE[0] * self.TE[1],
                                self.vals(self.TE_ref[0]) * self.vals(self.TE_ref[1]))
        self.assertEvalAllEqual(self.TE[0] * 0.5, self.vals(self.TE_ref[0]) * 0.5)
        self.assertEvalAllEqual(0.5 * self.TE[1], 0.5 * self.vals(self.TE_ref[1]))
        self.assertEvalAllEqual(self.TE[0] * 0.5j, self.vals(self.TE_ref[0]) * 0.5j)
        self.assertEvalAllEqual(0.5j * self.TE[1], 0.5j * self.vals(self.TE_ref[1]))

    def test_division(self):
        self.assertEvalAllEqual(self.TE[4] / 0.5, self.vals(self.TE_ref[4]) / 0.5)
        self.assertEvalAllEqual(self.TE[4] / 0.5j, self.vals(self.TE_ref[4]) / 0.5j)

if __name__ == '__main__':
    unittest.main()
