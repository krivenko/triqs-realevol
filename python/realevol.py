from realevol_core import texpr_r, texpr_c
from realevol_core import rmesh
from realevol_core import is_constant, is_zero

class texpr(object):
    """Factory object for both real and complex time-dependent expressions"""
    def __new__(texpr,*args):
        arity = {0: texpr_r(),
                 1: texpr_c(args[0]) if isinstance(args[0],complex) else texpr_r(args[0]),
                 2: texpr_c(*args)}
        if len(args) not in arity: raise TypeError("Too many arguments to construct texpr")
        return arity[len(args)]