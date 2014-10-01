from realevol_core import texpr_r, texpr_c
from realevol_core import rmesh
from realevol_core import Operator_r, c_r, c_dag_r, n_r
from realevol_core import Operator_c, c_c, c_dag_c, n_c
from realevol_core import is_constant, is_zero, dagger

class texpr(object):
    """Factory object for both real and complex time-dependent expressions"""
    def __new__(texpr,*args):
        arity = {0: texpr_r(),
                 1: texpr_c(args[0]) if isinstance(args[0],complex) else texpr_r(args[0]),
                 2: texpr_c(*args)}
        if len(args) not in arity: raise TypeError("Too many arguments to construct texpr")
        return arity[len(args)]

_complex_ops = False
def use_complex_operators(use):
    """Switch between real- and complex-valued fermionic operators"""
    global _complex_ops
    _complex_ops = use

class Operator(object):
    """Factory object for both real and complex fermionic operators"""
    def __new__(Operator):
        global _complex_ops
        return Operator_c() if _complex_ops else Operator_r()

def c(ind):
    global _complex_ops
    return c_c(ind) if _complex_ops else c_r(ind)
def c_dag(ind):
    global _complex_ops
    return c_dag_c(ind) if _complex_ops else c_dag_r(ind)
def n(ind):
    global _complex_ops
    return n_c(ind) if _complex_ops else n_r(ind)