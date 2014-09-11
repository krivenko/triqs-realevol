from pytriqs.applications.realevol.texpr import *
from math import sin, sqrt, pi


te1 = texpr("t^2")
te2 = texpr("t + sin(_pi/2)")
te3 = texpr("sqrt(9.0) + 1.5")

print "Check whether the expressions really depend on time"
print is_constant(te1)
print is_constant(te2)
print is_constant(te3)
if is_constant(te1) or is_constant(te2) or not is_constant(te3): exit(-1)

print "Check correctness of numerical expressions"
def check(x,y):
    print x,y
    if abs(x-y) >= 1e-10: exit(-1)

for t in [0, 0.1, 10, 55]:
    print "========","t =",t,"========"
    te1_ref = t**2
    te2_ref = t+sin(pi/2)
    te3_ref = sqrt(9.0)+1.5
    check(te1(t),te1_ref)
    check(te2(t),te2_ref)
    check(te3(t),te3_ref)

    print "Unary minus"
    mte2 = -te2
    check(mte2(t),-te2_ref) 

    print "Addition of expressions"
    te1pte2 = te1 + te2
    te1phalf = te1 + 0.5
    halfpte2 = 0.5 + te2
    check(te1pte2(t),te1_ref+te2_ref)
    check(te1phalf(t),te1_ref+0.5)
    check(halfpte2(t),0.5+te2_ref)

    print "Subtraction of expressions"
    te1mte2 = te1 - te2
    te1mhalf = te1 - 0.5
    halfmte2 = 0.5 - te2
    check(te1mte2(t),te1_ref-te2_ref)
    check(te1mhalf(t),te1_ref-0.5)
    check(halfmte2(t),0.5-te2_ref)

    print "Multiplication of expressions"
    te1ppte2 = te1 * te2
    te1pphalf = te1 * 0.5
    halfppte2 = 0.5 * te2
    check(te1ppte2(t),te1_ref*te2_ref)
    check(te1pphalf(t),te1_ref*0.5)
    check(halfppte2(t),0.5*te2_ref)

    print "Division of expressions by a scalar"
    te1dhalf = te1 / 0.5
    check(te1dhalf(t),te1_ref/0.5)