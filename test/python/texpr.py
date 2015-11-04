from pytriqs.applications.impurity_solvers.realevol.texpr import texpr, is_constant, is_zero, conj
from math import sin, sqrt, pi
from numpy import array

TE = [texpr("t^2"),
      texpr("t + sin(pi/2)"),
      texpr("sqrt(9.0) + 1.5"),
      texpr("t^2",1.0),
      texpr("t + sin(pi/2)","t^3"),
      texpr("sqrt(9.0) + 1.5","2.9"),
      texpr(1.9+2.8j)]

TE_ref = [lambda t: t**2,
          lambda t: t+sin(pi/2),
          lambda t: sqrt(9.0)+1.5,
          lambda t: t**2+1j,
          lambda t: t+sin(pi/2)+((t**3)*1j),
          lambda t: sqrt(9.0)+1.5+2.9j,
          lambda t: 1.9+2.8j]

def tab_print(title, A):
  print ("%20s "*(A.shape[0] + 1)) % ((title,) + tuple(a for a in A))

print "Check whether the expressions really depend on time"
for n,te in enumerate(TE): print "te[%i]:" % n, is_constant(te)

print "Check correctness of numerical expressions"

T = [0, 0.25, 16, 64]
tab_print("t:", array(T))

vals = lambda f: array(map(f,T),dtype=complex)

print "Evaluation"
for n, (te, te_ref) in enumerate(zip(TE,TE_ref)):
    tab_print("te[%i]" % n, vals(te))
    tab_print("te_ref[%i]" % n, vals(te_ref))

print "Unary minus"
tab_print("-te[1]", vals(-TE[1]))
tab_print("-te_ref[1]", -vals(TE_ref[1]))
tab_print("-te[4]", vals(-TE[4]))
tab_print("-te_ref[4]", -vals(TE_ref[4]))

print "Addition of expressions"
tab_print("te[0]+te[1]", vals(TE[0]) + vals(TE[1]))
tab_print("te_ref[0]+te_ref[1]", vals(TE_ref[0]) + vals(TE_ref[1]))
tab_print("te[0]+0.5", vals(TE[0] + 0.5))
tab_print("te_ref[0]+0.5", vals(TE_ref[0]) + 0.5)
tab_print("0.5+te[1]", vals(0.5 + TE[1]))
tab_print("0.5+te_ref[1]", 0.5 + vals(TE_ref[1]))
tab_print("te[0]+0.5j", vals(TE[0] + 0.5j))
tab_print("te_ref[0]+0.5j", vals(TE_ref[0]) + 0.5j)
tab_print("0.5j+te[1]", vals(0.5j + TE[1]))
tab_print("0.5j+te_ref[1]", 0.5j + vals(TE_ref[1]))

print "Subtraction of expressions"
tab_print("te[0]-te[1]", vals(TE[0]) - vals(TE[1]))
tab_print("te_ref[0]-te_ref[1]", vals(TE_ref[0]) - vals(TE_ref[1]))
tab_print("te[0]-0.5", vals(TE[0] - 0.5))
tab_print("te_ref[0]-0.5", vals(TE_ref[0]) - 0.5)
tab_print("0.5-te[1]", vals(0.5 - TE[1]))
tab_print("0.5-te_ref[1]", 0.5 - vals(TE_ref[1]))
tab_print("te[0]-0.5j", vals(TE[0] - 0.5j))
tab_print("te_ref[0]-0.5j", vals(TE_ref[0]) - 0.5j)
tab_print("0.5j-te[1]", vals(0.5j - TE[1]))
tab_print("0.5j-te_ref[1]", 0.5j - vals(TE_ref[1]))

print "Multiplication of expressions"
tab_print("te[0]*te[1]", vals(TE[0]) * vals(TE[1]))
tab_print("te_ref[0]*te_ref[1]", vals(TE_ref[0]) * vals(TE_ref[1]))
tab_print("te[0]*0.5", vals(TE[0] * 0.5))
tab_print("te_ref[0]*0.5", vals(TE_ref[0]) * 0.5)
tab_print("0.5*te[1]", vals(0.5 * TE[1]))
tab_print("0.5*te_ref[1]", 0.5 * vals(TE_ref[1]))
tab_print("te[0]*0.5j", vals(TE[0] * 0.5j))
tab_print("te_ref[0]*0.5j", vals(TE_ref[0]) * 0.5j)
tab_print("0.5j*te[1]", vals(0.5j * TE[1]))
tab_print("0.5j*te_ref[1]", 0.5j * vals(TE_ref[1]))

print "Division of expressions by a scalar"
tab_print("te[4]/0.5", vals(TE[4]/0.5))
tab_print("te_ref[4]/0.5", vals(TE_ref[4])/0.5)
tab_print("te[4]/0.5j", vals(TE[4]/0.5j))
tab_print("te_ref[4]/0.5j", vals(TE_ref[4])/0.5j)
