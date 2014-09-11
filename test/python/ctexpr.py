from pytriqs.applications.realevol.ctexpr import *
from pytriqs.applications.realevol.texpr import texpr as rtexpr
from cmath import sin, sqrt, pi

te1 = texpr("t^2",1.0)
te2 = texpr("t + sin(_pi/2)","t^3")
te3 = texpr("sqrt(9.0) + 1.5","2.9")
te4 = texpr(1.9+2.8j);

print "Check whether the expressions really depend on time"
print is_constant(te1)
print is_constant(te2)
print is_constant(te3)
print is_constant(te4)
if is_constant(te1) or is_constant(te2) or not is_constant(te3) or not is_constant(te4): exit(-1)

print "Check correctness of numerical expressions"
def check(x,y):
    print x,y
    if abs(x-y) >= 1e-10: exit(-1)

for t in [0, 0.1, 10, 55]:
    print "========","t =",t,"========"
    te1_ref = t**2+1j
    te2_ref = t+sin(pi/2)+((t**3)*1j)
    te3_ref = sqrt(9.0)+1.5+2.9j
    te4_ref = 1.9+2.8j
    check(te1(t),te1_ref)
    check(te2(t),te2_ref)
    check(te3(t),te3_ref)
    check(te4(t),te4_ref)

    print "Unary minus"
    mte2 = -te2
    check(mte2(t),-te2_ref) 

    print "Addition of expressions"
    e1pte2 = te1 + te2
    te1phalf = te1 + texpr(0.5)
    halfpte2 = texpr(0.5) + te2
    te1pihalf = te1 + texpr(0.5j)
    ihalfpte2 = texpr(0.5j) + te2
    check(e1pte2(t),te1_ref+te2_ref)
    check(te1phalf(t),te1_ref+0.5)
    check(halfpte2(t),0.5+te2_ref)
    check(te1pihalf(t),te1_ref+0.5j)
    check(ihalfpte2(t),0.5j+te2_ref)

    print "Subtraction of expressions"
    e1mte2 = te1 - te2
    te1mhalf = te1 - texpr(0.5)
    halfmte2 = texpr(0.5) - te2
    te1mihalf = te1 - texpr(0.5j)
    ihalfmte2 = texpr(0.5j) - te2
    check(e1mte2(t),te1_ref-te2_ref)
    check(te1mhalf(t),te1_ref-0.5)
    check(halfmte2(t),0.5-te2_ref)
    check(te1mihalf(t),te1_ref-0.5j)
    check(ihalfmte2(t),0.5j-te2_ref)

    print "Multiplication of expressions"
    te1ppte2 = te1 * te2
    te1pphalf = te1 * texpr(0.5)
    halfppte2 = texpr(0.5) * te2
    te1ppihalf = te1 * texpr(0.5j)
    ihalfppte2 = texpr(0.5j) * te2
    check(te1ppte2(t),te1_ref*te2_ref)
    check(te1pphalf(t),te1_ref*0.5)
    check(halfppte2(t),0.5*te2_ref)
    check(te1ppihalf(t),te1_ref*0.5j)
    check(ihalfppte2(t),0.5j*te2_ref)

    print "Division of expressions"
    te1dhalf = te1 / rtexpr(0.5)
    check(te1dhalf(t),te1_ref/0.5)