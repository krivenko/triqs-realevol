# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2021, I. Krivenko, M. Danilov, P. Kubiczek
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

from realevol.texpr import TExpr as te
from realevol.operators_texpr import *
from itertools import product

class test_operators_texpr(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.C = te("t^2") * c("",0)
        cls.Cd = c_dag("",1)
        cls.A = a("",0)
        cls.Ad = te("sin(t)") * a_dag("",1)

    def assertStrEqual(self, obj, s):
        self.assertEqual(str(obj), s)

    def test_no_indices(self):
        self.assertStrEqual(c() + c_dag() - n(), "1*C^+() + 1*C() + -1*C^+()C()")

    def test_commutators(self):
        for i,j in product(range(3), range(3)):
            self.assertStrEqual(c_dag("x",i)*c("x",j) + c("x",j)*c_dag("x",i),
                                "1" if i==j else "0")
            self.assertStrEqual(c_dag("x",i)*c("x",j) - c("x",j)*c_dag("x",i),
                                "%s2*C^+(x,%i)C(x,%i)" % ("-1 + " if i==j else "",i,j))
            self.assertStrEqual(a("x",i)*a_dag("x",j) - a_dag("x",j)*a("x",i),
                                "1" if i==j else "0")
            self.assertStrEqual(a("x",i)*a_dag("x",j) + a_dag("x",j)*a("x",i),
                                "%s2*A^+(x,%i)A(x,%i)" % ("1 + " if i==j else "",j,i))

    def test_canonical_ops(self):
        self.assertStrEqual(self.C, "t^2*C(,0)")
        self.assertStrEqual(self.Cd, "1*C^+(,1)")
        self.assertStrEqual(self.A, "1*A(,0)")
        self.assertStrEqual(self.Ad, "sin(t)*A^+(,1)")

    def test_minus(self):
        self.assertStrEqual(-self.C, "-(t^2)*C(,0)")
        self.assertStrEqual(-self.Cd, "-1*C^+(,1)")
        self.assertStrEqual(-self.A, "-1*A(,0)")
        self.assertStrEqual(-self.Ad, "-(sin(t))*A^+(,1)")

    def test_addition(self):
        self.assertStrEqual(self.C + 2.0, "2 + t^2*C(,0)")
        self.assertStrEqual(self.Cd + 2.0, "2 + 1*C^+(,1)")
        self.assertStrEqual(self.A + 2.0, "2 + 1*A(,0)")
        self.assertStrEqual(self.Ad + 2.0, "2 + sin(t)*A^+(,1)")
        self.assertStrEqual(self.C + 2.0j, "(0,2) + t^2*C(,0)")
        self.assertStrEqual(self.Cd + 2.0j, "(0,2) + 1*C^+(,1)")
        self.assertStrEqual(self.A + 2.0j, "(0,2) + 1*A(,0)")
        self.assertStrEqual(self.Ad + 2.0j, "(0,2) + sin(t)*A^+(,1)")
        self.assertStrEqual(self.C + te("2*t"), "2*t + t^2*C(,0)")
        self.assertStrEqual(self.Cd + te("2*t"), "2*t + 1*C^+(,1)")
        self.assertStrEqual(self.A + te("2*t"), "2*t + 1*A(,0)")
        self.assertStrEqual(self.Ad + te("2*t"), "2*t + sin(t)*A^+(,1)")
        self.assertStrEqual(self.C + 1j*te("2*t"), "(0,2*t) + t^2*C(,0)")
        self.assertStrEqual(self.Cd + 1j*te("2*t"), "(0,2*t) + 1*C^+(,1)")
        self.assertStrEqual(self.A + 1j*te("2*t"), "(0,2*t) + 1*A(,0)")
        self.assertStrEqual(self.Ad + 1j*te("2*t"), "(0,2*t) + sin(t)*A^+(,1)")
        self.assertStrEqual(2.0 + self.C, "2 + t^2*C(,0)")
        self.assertStrEqual(2.0 + self.Cd, "2 + 1*C^+(,1)")
        self.assertStrEqual(2.0 + self.A, "2 + 1*A(,0)")
        self.assertStrEqual(2.0 + self.Ad, "2 + sin(t)*A^+(,1)")
        self.assertStrEqual(2.0j + self.C, "(0,2) + t^2*C(,0)")
        self.assertStrEqual(2.0j + self.Cd, "(0,2) + 1*C^+(,1)")
        self.assertStrEqual(2.0j + self.A, "(0,2) + 1*A(,0)")
        self.assertStrEqual(2.0j + self.Ad, "(0,2) + sin(t)*A^+(,1)")
        self.assertStrEqual(te("2*t") + self.C, "2*t + t^2*C(,0)")
        self.assertStrEqual(te("2*t") + self.Cd, "2*t + 1*C^+(,1)")
        self.assertStrEqual(te("2*t") + self.A, "2*t + 1*A(,0)")
        self.assertStrEqual(te("2*t") + self.Ad, "2*t + sin(t)*A^+(,1)")
        self.assertStrEqual(1j*te("2*t") + self.C, "(0,2*t) + t^2*C(,0)")
        self.assertStrEqual(1j*te("2*t") + self.Cd, "(0,2*t) + 1*C^+(,1)")
        self.assertStrEqual(1j*te("2*t") + self.A, "(0,2*t) + 1*A(,0)")
        self.assertStrEqual(1j*te("2*t") + self.Ad, "(0,2*t) + sin(t)*A^+(,1)")

        self.assertStrEqual(self.C + self.Cd + self.A + self.Ad,
                            "1*C^+(,1) + t^2*C(,0) + sin(t)*A^+(,1) + 1*A(,0)")

    def test_subtraction(self):
        self.assertStrEqual(self.C - 2.0, "-2 + t^2*C(,0)")
        self.assertStrEqual(self.Cd - 2.0, "-2 + 1*C^+(,1)")
        self.assertStrEqual(self.A - 2.0, "-2 + 1*A(,0)")
        self.assertStrEqual(self.Ad - 2.0, "-2 + sin(t)*A^+(,1)")
        self.assertStrEqual(self.C - 2.0j, "(0,-2) + t^2*C(,0)")
        self.assertStrEqual(self.Cd - 2.0j, "(0,-2) + 1*C^+(,1)")
        self.assertStrEqual(self.A - 2.0j, "(0,-2) + 1*A(,0)")
        self.assertStrEqual(self.Ad - 2.0j, "(0,-2) + sin(t)*A^+(,1)")
        self.assertStrEqual(self.C - te("2*t"), "-(2*t) + t^2*C(,0)")
        self.assertStrEqual(self.Cd - te("2*t"), "-(2*t) + 1*C^+(,1)")
        self.assertStrEqual(self.A - te("2*t"), "-(2*t) + 1*A(,0)")
        self.assertStrEqual(self.Ad - te("2*t"), "-(2*t) + sin(t)*A^+(,1)")
        self.assertStrEqual(self.C - 1j*te("2*t"), "(0,-(2*t)) + t^2*C(,0)")
        self.assertStrEqual(self.Cd - 1j*te("2*t"), "(0,-(2*t)) + 1*C^+(,1)")
        self.assertStrEqual(self.A - 1j*te("2*t"), "(0,-(2*t)) + 1*A(,0)")
        self.assertStrEqual(self.Ad - 1j*te("2*t"), "(0,-(2*t)) + sin(t)*A^+(,1)")
        self.assertStrEqual(2.0 - self.C, "2 + -(t^2)*C(,0)")
        self.assertStrEqual(2.0 - self.Cd, "2 + -1*C^+(,1)")
        self.assertStrEqual(2.0 - self.A, "2 + -1*A(,0)")
        self.assertStrEqual(2.0 - self.Ad, "2 + -(sin(t))*A^+(,1)")
        self.assertStrEqual(2.0j - self.C, "(0,2) + -(t^2)*C(,0)")
        self.assertStrEqual(2.0j - self.Cd, "(0,2) + -1*C^+(,1)")
        self.assertStrEqual(2.0j - self.A, "(0,2) + -1*A(,0)")
        self.assertStrEqual(2.0j - self.Ad, "(0,2) + -(sin(t))*A^+(,1)")
        self.assertStrEqual(te("2*t") - self.C, "2*t + -(t^2)*C(,0)")
        self.assertStrEqual(te("2*t") - self.Cd, "2*t + -1*C^+(,1)")
        self.assertStrEqual(te("2*t") - self.A, "2*t + -1*A(,0)")
        self.assertStrEqual(te("2*t") - self.Ad, "2*t + -(sin(t))*A^+(,1)")
        self.assertStrEqual(1j*te("2*t") - self.C, "(0,2*t) + -(t^2)*C(,0)")
        self.assertStrEqual(1j*te("2*t") - self.Cd, "(0,2*t) + -1*C^+(,1)")
        self.assertStrEqual(1j*te("2*t") - self.A, "(0,2*t) + -1*A(,0)")
        self.assertStrEqual(1j*te("2*t") - self.Ad, "(0,2*t) + -(sin(t))*A^+(,1)")

        self.assertStrEqual(self.C - self.Cd - self.A - self.Ad,
                            "-1*C^+(,1) + t^2*C(,0) + -(sin(t))*A^+(,1) + -1*A(,0)")

    def test_multiplication(self):
        self.assertStrEqual(self.C * 3.0, "(t^2)*(3)*C(,0)")
        self.assertStrEqual(self.Cd * 3.0, "3*C^+(,1)")
        self.assertStrEqual(self.A * 3.0, "3*A(,0)")
        self.assertStrEqual(self.Ad * 3.0, "(sin(t))*(3)*A^+(,1)")
        self.assertStrEqual(self.C * 3.0j, "(0,(t^2)*(3))*C(,0)")
        self.assertStrEqual(self.Cd * 3.0j, "(0,3)*C^+(,1)")
        self.assertStrEqual(self.A * 3.0j, "(0,3)*A(,0)")
        self.assertStrEqual(self.Ad * 3.0j, "(0,(sin(t))*(3))*A^+(,1)")
        self.assertStrEqual(self.C * te("2*t"), "(t^2)*(2*t)*C(,0)")
        self.assertStrEqual(self.Cd * te("2*t"), "2*t*C^+(,1)")
        self.assertStrEqual(self.A * te("2*t"), "2*t*A(,0)")
        self.assertStrEqual(self.Ad * te("2*t"), "(sin(t))*(2*t)*A^+(,1)")
        self.assertStrEqual(self.C * 1j*te("2*t"), "(0,(t^2)*(2*t))*C(,0)")
        self.assertStrEqual(self.Cd * 1j*te("2*t"), "(0,2*t)*C^+(,1)")
        self.assertStrEqual(self.A * 1j*te("2*t"), "(0,2*t)*A(,0)")
        self.assertStrEqual(self.Ad * 1j*te("2*t"), "(0,(sin(t))*(2*t))*A^+(,1)")
        self.assertStrEqual(3.0 * self.C, "(t^2)*(3)*C(,0)")
        self.assertStrEqual(3.0 * self.Cd, "3*C^+(,1)")
        self.assertStrEqual(3.0 * self.A, "3*A(,0)")
        self.assertStrEqual(3.0 * self.Ad, "(sin(t))*(3)*A^+(,1)")
        self.assertStrEqual(3.0j * self.C, "(0,(t^2)*(3))*C(,0)")
        self.assertStrEqual(3.0j * self.Cd, "(0,3)*C^+(,1)")
        self.assertStrEqual(3.0j * self.A, "(0,3)*A(,0)")
        self.assertStrEqual(3.0j * self.Ad, "(0,(sin(t))*(3))*A^+(,1)")
        self.assertStrEqual(te("2*t") * self.C, "(t^2)*(2*t)*C(,0)")
        self.assertStrEqual(te("2*t") * self.Cd, "2*t*C^+(,1)")
        self.assertStrEqual(te("2*t") * self.A, "2*t*A(,0)")
        self.assertStrEqual(te("2*t") * self.Ad, "(sin(t))*(2*t)*A^+(,1)")
        self.assertStrEqual(1j*te("2*t") * self.C, "(0,(t^2)*(2*t))*C(,0)")
        self.assertStrEqual(1j*te("2*t") * self.Cd, "(0,2*t)*C^+(,1)")
        self.assertStrEqual(1j*te("2*t") * self.A, "(0,2*t)*A(,0)")
        self.assertStrEqual(1j*te("2*t") * self.Ad, "(0,(sin(t))*(2*t))*A^+(,1)")

        self.assertStrEqual((c("",2) + c_dag("",2) + a("",2) + a_dag("",2)) *
                            (c("",2) - c_dag("",2) + a("",2) - a_dag("",2)),
                            "-2 + 2*C^+(,2)C(,2) + -2*C^+(,2)A^+(,2) + 2*C(,2)A(,2) + "
                            "-1*[A^+(,2)]^2 + 1*[A(,2)]^2")

    def test_N3(self):
        expr = n("up",0) * a(0,0) + n("dn",0) * a_dag(0,0)
        self.assertStrEqual(expr, "1*C^+(dn,0)C(dn,0)A^+(0,0) + 1*C^+(up,0)C(up,0)A(0,0)")
        expr = expr * expr * expr
        self.assertStrEqual(expr,
                            "3*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)A^+(0,0) + "
                            "3*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)A(0,0) + "
                            "1*C^+(dn,0)C(dn,0)[A^+(0,0)]^3 + "
                            "1*C^+(up,0)C(up,0)[A(0,0)]^3 + "
                            "3*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)[A^+(0,0)]^2A(0,0) + "
                            "3*C^+(dn,0)C^+(up,0)C(up,0)C(dn,0)A^+(0,0)[A(0,0)]^2")

    def test_dagger(self):
        X = te("t^2","sin(t)")*c_dag("",1) * c_dag("",2) * c("",3) * c("",4) * a_dag("",5) * a("",6)
        self.assertStrEqual(X, "(-(t^2),-(sin(t)))*C^+(,1)C^+(,2)C(,4)C(,3)A^+(,5)A(,6)")
        self.assertStrEqual(dagger(X), "(-(t^2),-(-(sin(t))))*C^+(,3)C^+(,4)C(,2)C(,1)A^+(,6)A(,5)")

if __name__ == '__main__':
    unittest.main()
