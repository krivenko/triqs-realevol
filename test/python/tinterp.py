from realevol.tinterp import TInterp, is_constant, is_zero, conj
from pytriqs.gf.local import MeshReTime
from numpy import array, where

m = MeshReTime(0, 10, 11)

TI = [TInterp(),
      TInterp(m, array([(t*t).real for t in m])),
      TInterp(m, array([(t + 1).real for t in m])),
      TInterp(m, 4.5),
      TInterp(m, array([(t*t).real for t in m]), 1.0),
      TInterp(m, array([(t + 1).real for t in m]), array([(t*t).real for t in m])),
      TInterp(m, 4.5, 2.9),
      TInterp(m, 1.9+2.8j),
      TInterp(m, array([t*t + (t + 1)*1j for t in m]))]

T = array([0, 0.1, 1, 5.5, 10])
TI_res = [array([0, 0, 0, 0, 0]),
          array([0, 0.1, 1, 30.5, 100]),
          array([1, 1.1, 2, 6.5, 11]),
          array([4.5, 4.5, 4.5, 4.5, 4.5]),
          array([1.0j, 0.1+1.0j, 1+1.0j, 30.5+1.0j, 100+1.0j]),
          array([1, 1.1+0.1j, 2+1j, 6.5+30.5j, 11+100j]),
          array([4.5+2.9j, 4.5+2.9j, 4.5+2.9j, 4.5+2.9j, 4.5+2.9j]),
          array([1.9+2.8j, 1.9+2.8j, 1.9+2.8j, 1.9+2.8j, 1.9+2.8j]),
          array([1j, 0.1+1.1j, 1+2j, 30.5+6.5j, 100+11j])]

# Check constantness of interpolators
assert map(is_constant, TI) == [True, False, False, True, False, False, True, True, False]
assert is_zero(TI[0])
assert not is_zero(TI[1])

# Check conj()
assert conj(TI[3]) == TInterp(m, 4.5)
assert conj(TI[5]) == TInterp(m, array([t + 1 - t*t*1j for t in m]))

def check_eval(ti, ti_res):
    return all(abs(array(map(ti, T)) - ti_res) < 1e-12)

# Check correctness of evaluated values
for i in range(len(TI)):
  assert check_eval(TI[i], TI_res[i])

# Unary minus
assert check_eval(-TI[2], -TI_res[2])
assert check_eval(-TI[5], -TI_res[5])

# Addition
assert check_eval(TI[1] + TI[2], TI_res[1] + TI_res[2])
assert check_eval(TI[1] + 0.5, TI_res[1] + 0.5)
assert check_eval(0.5 + TI[2], 0.5 + TI_res[2])
assert check_eval(TI[1] + 0.5j, TI_res[1] + 0.5j)
assert check_eval(0.5j + TI[2], 0.5j + TI_res[2])

# Subtraction
assert check_eval(TI[1] - TI[2], TI_res[1] - TI_res[2])
assert check_eval(TI[1] - 0.5, TI_res[1] - 0.5)
assert check_eval(0.5 - TI[2], 0.5 - TI_res[2])
assert check_eval(TI[1] - 0.5j, TI_res[1] - 0.5j)
assert check_eval(0.5j - TI[2], 0.5j - TI_res[2])

# Multiplication
assert check_eval(TI[1] * TI[2], array([0, 0.2, 2, 201, 1100]))
assert check_eval(TI[5] * TI[6], array([4.5+2.9j, 4.66+3.64j, 6.1+10.3j, -59.2+156.1j, -240.5+481.9j]))
assert check_eval(TI[1] * 0.5, TI_res[1] * 0.5)
assert check_eval(0.5 * TI[2], 0.5 * TI_res[2])
assert check_eval(TI[1] * 0.5j, TI_res[1] * 0.5j)
assert check_eval(0.5j * TI[2], 0.5j * TI_res[2])

# Division
assert check_eval(TI[1] / 0.5, TI_res[1] / 0.5)
assert check_eval(TI[1] / 0.5j, TI_res[1] / 0.5j)
