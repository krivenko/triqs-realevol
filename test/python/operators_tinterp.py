# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2025, I. Krivenko, M. Danilov, P. Kubiczek
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

from realevol.tinterp import TInterp as ti
from realevol.operators_tinterp import *
from triqs.gf import MeshReTime
from itertools import product
from numpy import array

class test_operators_tinterp(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.m = MeshReTime(0, 1, 6)

        cls.ti1 = ti(cls.m, array([.0, 0.2, 0.4, 0.6, 0.8, 1.0]))
        cls.ti2 = ti(cls.m, array([.0, 0.1, 0.3, 0.5, 0.7, 0.9]))
        cls.ti3 = ti(cls.m, array([.0, 0.2, 0.4, 0.6, 0.4, 0.2]))

        cls.C = cls.ti1 * c("",0)
        cls.Cd = c_dag("",1);
        cls.A = a("",0)
        cls.Ad = cls.ti2 * a_dag("",1)

    def assertStrEqual(self, obj, s):
        self.assertEqual(str(obj), s)

    def test_no_indices(self):
        self.assertStrEqual(c() + c_dag() - n(),
                            "(1,0)*C^+() + (1,0)*C() + (-1,-0)*C^+()C()")

    def test_commutators(self):
        for i,j in product(range(3), range(3)):
            self.assertStrEqual(c_dag("x",i)*c("x",j) + c("x",j)*c_dag("x",i),
                                "(1,0)" if i==j else "0")
            self.assertStrEqual(c_dag("x",i)*c("x",j) - c("x",j)*c_dag("x",i),
                                "%s(2,0)*C^+(x,%i)C(x,%i)" % ("(-1,-0) + " if i==j else "",i,j))
            self.assertStrEqual(a("x",i)*a_dag("x",j) - a_dag("x",j)*a("x",i),
                                "(1,0)" if i==j else "0")
            self.assertStrEqual(a("x",i)*a_dag("x",j) + a_dag("x",j)*a("x",i),
                                "%s(2,0)*A^+(x,%i)A(x,%i)" % ("(1,0) + " if i==j else "",j,i))

    def test_canonical_ops(self):
        self.assertStrEqual(self.C, "ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd, "(1,0)*C^+(,1)")
        self.assertStrEqual(self.A, "(1,0)*A(,0)")
        self.assertStrEqual(self.Ad, "ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")

    def test_minus(self):
        self.assertStrEqual(-self.C, "ti([0,1]->[(-0,-0),...,(-1,-0)])*C(,0)")
        self.assertStrEqual(-self.Cd, "(-1,-0)*C^+(,1)")
        self.assertStrEqual(-self.A, "(-1,-0)*A(,0)")
        self.assertStrEqual(-self.Ad, "ti([0,1]->[(-0,-0),...,(-0.9,-0)])*A^+(,1)")

    def test_addition(self):
        self.assertStrEqual(self.C + 2.0, "(2,0) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd + 2.0, "(2,0) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A + 2.0, "(2,0) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad + 2.0, "(2,0) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(self.C + 2.0j, "(0,2) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd + 2.0j, "(0,2) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A + 2.0j, "(0,2) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad + 2.0j, "(0,2) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(self.C + self.ti3,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd + self.ti3,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A + self.ti3,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad + self.ti3,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(self.C + 1j*self.ti3,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd + 1j*self.ti3,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A + 1j*self.ti3,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad + 1j*self.ti3,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(2.0 + self.C, "(2,0) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(2.0 + self.Cd, "(2,0) + (1,0)*C^+(,1)")
        self.assertStrEqual(2.0 + self.A, "(2,0) + (1,0)*A(,0)")
        self.assertStrEqual(2.0 + self.Ad, "(2,0) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(2.0j + self.C, "(0,2) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(2.0j + self.Cd, "(0,2) + (1,0)*C^+(,1)")
        self.assertStrEqual(2.0j + self.A, "(0,2) + (1,0)*A(,0)")
        self.assertStrEqual(2.0j + self.Ad, "(0,2) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(self.ti3 + self.C,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.ti3 + self.Cd,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.ti3 + self.A,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + (1,0)*A(,0)")
        self.assertStrEqual(self.ti3 + self.Ad,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(1j*self.ti3 + self.C,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(1j*self.ti3 + self.Cd,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + (1,0)*C^+(,1)")
        self.assertStrEqual(1j*self.ti3 + self.A,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + (1,0)*A(,0)")
        self.assertStrEqual(1j*self.ti3 + self.Ad,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")

        self.assertStrEqual(self.C + self.Cd + self.A + self.Ad,
                            "(1,0)*C^+(,1) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)"
                            " + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1) + (1,0)*A(,0)")

    def test_subtraction(self):
        self.assertStrEqual(self.C - 2.0, "(-2,-0) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd - 2.0, "(-2,-0) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A - 2.0, "(-2,-0) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad - 2.0, "(-2,-0) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(self.C - 2.0j, "(-0,-2) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd - 2.0j, "(-0,-2) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A - 2.0j, "(-0,-2) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad - 2.0j, "(-0,-2) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(self.C - self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0.2,-0)]) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)",)
        self.assertStrEqual(self.Cd - self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0.2,-0)]) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A - self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0.2,-0)]) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad - self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0.2,-0)]) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(self.C - 1j*self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0,-0.2)]) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)")
        self.assertStrEqual(self.Cd - 1j*self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0,-0.2)]) + (1,0)*C^+(,1)")
        self.assertStrEqual(self.A - 1j*self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0,-0.2)]) + (1,0)*A(,0)")
        self.assertStrEqual(self.Ad - 1j*self.ti3,
                            "ti([0,1]->[(-0,-0),...,(-0,-0.2)]) + ti([0,1]->[(0,0),...,(0.9,0)])*A^+(,1)")
        self.assertStrEqual(2.0 - self.C, "(2,0) + ti([0,1]->[(-0,-0),...,(-1,-0)])*C(,0)")
        self.assertStrEqual(2.0 - self.Cd, "(2,0) + (-1,-0)*C^+(,1)")
        self.assertStrEqual(2.0 - self.A, "(2,0) + (-1,-0)*A(,0)")
        self.assertStrEqual(2.0 - self.Ad, "(2,0) + ti([0,1]->[(-0,-0),...,(-0.9,-0)])*A^+(,1)")
        self.assertStrEqual(2.0j - self.C, "(0,2) + ti([0,1]->[(-0,-0),...,(-1,-0)])*C(,0)")
        self.assertStrEqual(2.0j - self.Cd, "(0,2) + (-1,-0)*C^+(,1)")
        self.assertStrEqual(2.0j - self.A, "(0,2) + (-1,-0)*A(,0)")
        self.assertStrEqual(2.0j - self.Ad, "(0,2) + ti([0,1]->[(-0,-0),...,(-0.9,-0)])*A^+(,1)")
        self.assertStrEqual(self.ti3 - self.C,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + ti([0,1]->[(-0,-0),...,(-1,-0)])*C(,0)")
        self.assertStrEqual(self.ti3 - self.Cd,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + (-1,-0)*C^+(,1)")
        self.assertStrEqual(self.ti3 - self.A,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + (-1,-0)*A(,0)")
        self.assertStrEqual(self.ti3 - self.Ad,
                            "ti([0,1]->[(0,0),...,(0.2,0)]) + ti([0,1]->[(-0,-0),...,(-0.9,-0)])*A^+(,1)")
        self.assertStrEqual(1j*self.ti3 - self.C,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + ti([0,1]->[(-0,-0),...,(-1,-0)])*C(,0)")
        self.assertStrEqual(1j*self.ti3 - self.Cd,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + (-1,-0)*C^+(,1)")
        self.assertStrEqual(1j*self.ti3 - self.A,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + (-1,-0)*A(,0)")
        self.assertStrEqual(1j*self.ti3 - self.Ad,
                            "ti([0,1]->[(0,0),...,(0,0.2)]) + ti([0,1]->[(-0,-0),...,(-0.9,-0)])*A^+(,1)")

        self.assertStrEqual(self.C - self.Cd - self.A - self.Ad,
                            "(-1,-0)*C^+(,1) + ti([0,1]->[(0,0),...,(1,0)])*C(,0)"
                            " + ti([0,1]->[(-0,-0),...,(-0.9,-0)])*A^+(,1) + (-1,-0)*A(,0)")

    def test_multiplication(self):
        self.assertStrEqual(self.C * 3.0, "ti([0,1]->[(0,0),...,(3,0)])*C(,0)")
        self.assertStrEqual(self.Cd * 3.0, "(3,0)*C^+(,1)")
        self.assertStrEqual(self.A * 3.0, "(3,0)*A(,0)")
        self.assertStrEqual(self.Ad * 3.0, "ti([0,1]->[(0,0),...,(2.7,0)])*A^+(,1)")
        self.assertStrEqual(self.C * 3.0j, "ti([0,1]->[(0,0),...,(0,3)])*C(,0)")
        self.assertStrEqual(self.Cd * 3.0j, "(0,3)*C^+(,1)")
        self.assertStrEqual(self.A * 3.0j, "(0,3)*A(,0)")
        self.assertStrEqual(self.Ad * 3.0j, "ti([0,1]->[(0,0),...,(0,2.7)])*A^+(,1)")
        self.assertStrEqual(self.C * self.ti3, "ti([0,1]->[(0,0),...,(0.2,0)])*C(,0)")
        self.assertStrEqual(self.Cd * self.ti3, "ti([0,1]->[(0,0),...,(0.2,0)])*C^+(,1)")
        self.assertStrEqual(self.A * self.ti3, "ti([0,1]->[(0,0),...,(0.2,0)])*A(,0)")
        self.assertStrEqual(self.Ad * self.ti3, "ti([0,1]->[(0,0),...,(0.18,0)])*A^+(,1)")
        self.assertStrEqual(self.C * 1j*self.ti3, "ti([0,1]->[(0,0),...,(0,0.2)])*C(,0)")
        self.assertStrEqual(self.Cd * 1j*self.ti3, "ti([0,1]->[(0,0),...,(0,0.2)])*C^+(,1)")
        self.assertStrEqual(self.A * 1j*self.ti3, "ti([0,1]->[(0,0),...,(0,0.2)])*A(,0)")
        self.assertStrEqual(self.Ad * 1j*self.ti3, "ti([0,1]->[(0,0),...,(0,0.18)])*A^+(,1)")
        self.assertStrEqual(3.0 * self.C, "ti([0,1]->[(0,0),...,(3,0)])*C(,0)")
        self.assertStrEqual(3.0 * self.Cd, "(3,0)*C^+(,1)")
        self.assertStrEqual(3.0 * self.A, "(3,0)*A(,0)")
        self.assertStrEqual(3.0 * self.Ad, "ti([0,1]->[(0,0),...,(2.7,0)])*A^+(,1)")
        self.assertStrEqual(3.0j * self.C, "ti([0,1]->[(0,0),...,(0,3)])*C(,0)")
        self.assertStrEqual(3.0j * self.Cd, "(0,3)*C^+(,1)")
        self.assertStrEqual(3.0j * self.A, "(0,3)*A(,0)")
        self.assertStrEqual(3.0j * self.Ad, "ti([0,1]->[(0,0),...,(0,2.7)])*A^+(,1)")
        self.assertStrEqual(self.ti3 * self.C, "ti([0,1]->[(0,0),...,(0.2,0)])*C(,0)")
        self.assertStrEqual(self.ti3 * self.Cd, "ti([0,1]->[(0,0),...,(0.2,0)])*C^+(,1)")
        self.assertStrEqual(self.ti3 * self.A, "ti([0,1]->[(0,0),...,(0.2,0)])*A(,0)")
        self.assertStrEqual(self.ti3 * self.Ad, "ti([0,1]->[(0,0),...,(0.18,0)])*A^+(,1)")
        self.assertStrEqual(1j*self.ti3 * self.C, "ti([0,1]->[(0,0),...,(0,0.2)])*C(,0)")
        self.assertStrEqual(1j*self.ti3 * self.Cd, "ti([0,1]->[(0,0),...,(0,0.2)])*C^+(,1)")
        self.assertStrEqual(1j*self.ti3 * self.A, "ti([0,1]->[(0,0),...,(0,0.2)])*A(,0)")
        self.assertStrEqual(1j*self.ti3 * self.Ad, "ti([0,1]->[(0,0),...,(0,0.18)])*A^+(,1)")

        self.assertStrEqual((c("",2) + c_dag("",2) + a("",2) + a_dag("",2)) *
                            (c("",2) - c_dag("",2) + a("",2) - a_dag("",2)),
                            "(-2,-0) + (2,0)*C^+(,2)C(,2) + (-2,-0)*C^+(,2)A^+(,2)"
                            " + (2,0)*C(,2)A(,2) + (-1,-0)*[A^+(,2)]^2 + (1,0)*[A(,2)]^2")

    def test_N3(self):
        expr = n("up",0) * a(0,0) + n("dn",0) * a_dag(0,0)
        self.assertStrEqual(expr, "(1,0)*C^+(dn,0)C(dn,0)A^+(0,0) + (1,0)*C^+(up,0)C(up,0)A(0,0)")
        expr = expr * expr * expr
        self.assertStrEqual(expr,
                            "(3,0)*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)A^+(0,0) + "
                            "(3,0)*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)A(0,0) + "
                            "(1,0)*C^+(dn,0)C(dn,0)[A^+(0,0)]^3 + "
                            "(1,0)*C^+(up,0)C(up,0)[A(0,0)]^3 + "
                            "(3,0)*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)[A^+(0,0)]^2A(0,0) + "
                            "(3,0)*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)A^+(0,0)[A(0,0)]^2")

    def test_dagger(self):
        X = (self.ti1 + 1j*self.ti2)*c_dag("",1) * c_dag("",2) * c("",3) * c("",4) * a_dag("",5) * a("",6)
        self.assertStrEqual(X, "ti([0,1]->[(0,0),...,(-1,-0.9)])*C^+(,1)C^+(,2)C(,4)C(,3)A^+(,5)A(,6)")
        self.assertStrEqual(dagger(X),
                            "ti([0,1]->[(0,-0),...,(-1,0.9)])*C^+(,3)C^+(,4)C(,2)C(,1)A^+(,6)A(,5)")

    def test_operator_stat(self):
        self.assertStrEqual(operator_stat(Operator()), "Boson")
        self.assertStrEqual(operator_stat(Operator(self.ti1)), "Boson")
        self.assertStrEqual(operator_stat(c_dag("",1)), "Fermion")
        self.assertStrEqual(operator_stat(c("",1)), "Fermion")
        self.assertStrEqual(operator_stat(a_dag("",1)), "Boson")
        self.assertStrEqual(operator_stat(a("",1)), "Boson")
        self.assertStrEqual(operator_stat(c_dag("",1) * c("",2) + c_dag("",2) * c("",1)), "Boson")
        self.assertStrEqual(operator_stat(c_dag("",1) * c("",2) * c("",3) + c_dag("",4)), "Fermion")
        self.assertStrEqual(operator_stat(a_dag("",1) * a("",2)), "Boson")
        self.assertStrEqual(operator_stat(a_dag("",1) * a("",2) * c("",3) + c_dag("",4)), "Fermion")
        with self.assertRaises(RuntimeError):
            operator_stat(a_dag("",1) * a("",2) * c("",3) + c_dag("",4) * c("",5))

if __name__ == '__main__':
    unittest.main()
