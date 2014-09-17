from pytriqs.applications.realevol.texpr import *
from pytriqs.applications.realevol.operators import *
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