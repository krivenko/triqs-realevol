from pytriqs.applications.realevol.realevol import texpr
from pytriqs.applications.realevol.realevol import use_complex_operators, Operator, c, c_dag, n, dagger
from itertools import product

def test_commutators(Cd,C):
    print
    print "Commutators:"
    for cdi, ci in product(Cd,C):
        print "[", cdi, ",", ci, "] =", cdi*ci - ci*cdi

def test_anticommutators(Cd,C):
    print
    print "Anticommutators:"
    for cdi, ci in product(Cd,C):
        print "{", cdi, ",", ci, "} =", cdi*ci + ci*cdi

print "Real-expression-valued operators"
print "====================================="

use_complex_operators(False)

# Commutation relations
C = [c(1), c(2), c(3)]
Cd = [c_dag(1), c_dag(2), c_dag(3)]

test_anticommutators(Cd,C)
test_commutators(Cd,C)

x = "t^2"*c(0)
y = c_dag(1)

print
print "Algebra:"
print "x =", x
print "y =", y

print "-x=", -x
print
print "x + 2.0 =", x + 2.0
print "2.0 + x =", 2.0 + x
print "x - 2.0 =", x - 2.0
print "2.0 - x =", 2.0 - x
print "3.0*y =", 3.0*y
print "y*3.0 =", y*3.0
print "x + 2.0*t =", x + "2*t"
print "2.0*t + x =", "2*t" + x
print "x - 2.0*t =", x - "2.0*t"
print "2.0*t - x =", "2.0*t" - x
print "3.0*t*y =", "3.0*t"*y
print "y*3.0*t =", y*"3.0*t"
print "x + y =", x + y
print "x - y =", x - y
print "(x + y)*(x - y) =", (x + y)*(x - y)

print
print "N^3:"
N = n("up") + n("dn")
N3 = N*N*N
print "N =", N
print "N^3 =", N3

X = c_dag(1) * c_dag(2) * c(3) * c(4)
print
print "X =", X
print "dagger(X) =", dagger(X)

print "Complex-expression-valued operators"
print "====================================="

use_complex_operators(True)

# Commutation relations
C = [c(1), c(2), c(3)]
Cd = [c_dag(1), c_dag(2), c_dag(3)]

test_anticommutators(Cd,C)
test_commutators(Cd,C)

x = "t^2"*c(0)
y = c_dag(1)

print
print "Algebra:"
print "x =", x
print "y =", y

print "-x=", -x
print
print "x + 2.0*I =", x + 2.0*1j
print "2.0*I + x =", 2.0*1j + x
print "x - 2.0*I =", x - 2.0*1j
print "2.0*I - x =", 2.0*1j - x
print "3.0*I*y =", (3.0*1j)*y
print "y*3.0*I =", y*(3.0*1j)
print "x + 2.0*t*I =", x + texpr("2*t",0)*1j
print "2.0*t*I + x =", texpr("2*t",0)*1j + x
print "x - 2.0*t*I =", x - texpr("2*t",0)*1j
print "2.0*t*I - x =", texpr("2.0*t",0)*1j - x
print "3.0*t*I*y =", texpr("3.0*t",0)*1j*y
print "y*3.0*I*t =", y*texpr("3.0*t",0)*1j
print "x + y =", x + y
print "x - y =", x - y
print "(x + y)*(x - y) =", (x + y)*(x - y)

print
print "N^3:"
N = n("up") + n("dn")
N3 = N*N*N
print "N =", N
print "N^3 =", N3

X = 1j*c_dag(1) * c_dag(2) * c(3) * c(4)
print
print "X =", X
print "dagger(X) =", dagger(X)