# realevol example: Hubbard atom coupled to one boson
# The coupling is slowly switched on at t = 0
#
# H(t) = eps0 * (n_up + n_dn) + U n_up n_dn + omega * a^+ a
#      + g(t) * (n_up + n_dn) * (a^+ + a),

from pytriqs.archive import HDFArchive  # HDF5 archive interface
from realevol.texpr import TExpr        # Time-dependent expressions
from realevol.operators_texpr import *  # Time-dependent second quantization operators
from realevol.gf_retime import *        # Green's function of two times
from realevol.init_state import *       # Routines to prepare initial states
from realevol.realevol import *         # Solver object
import pytriqs.utility.mpi as mpi       # MPI utilities

# Model parameters
T = 0.5                         # Temperature of the equilibrium initial state
U = 2.0                         # Hubbard repulsion
eps0 = -U/2                     # Position of local level
omega = 4.0                     # Energy of the bosonic excitation
g = TExpr("0.1*(1-exp(-5*t))")  # Time-dependent boson/fermion coupling constant

# Structure of greater/lesser Green's functions: two 1x1 blocks
# The blocks correspond to spin projections, which are not mixed by H(t)
gf_struct = {'dn' : [0], 'up' : [0]}

# Indices of susceptibility functions
# Solver will measure all <n_x n_y>, where x, y = (dn,0), (up,0) (4 components in total)
chi_indices = [('dn',0), ('up',0)]

# Maximum observation time
t_max = 1.0
# Number of slices in the time grid
n_t = 51

# Hamiltonian at t = 0
H0 = eps0 * (n('up',0) + n('dn',0)) + U*n('up',0)*n('dn',0)
H0 = H0 + omega * a_dag('B',0) * a('B',0)

# Set of all fermionic degrees of freedom
fermion_indices = set([('dn',0), ('up',0)])
# Set of all bosonic degrees of freedom
boson_indices = set([('B',0)])
# We allow up to 2^3 - 1 = 7 excitations in bosonic mode ('B',0)
bits_per_boson = {('B',0) : 3}

# Prepare equilibrium initial state
init_state = make_equilibrium_init_state(H0,
                                         temperature = T,
                                         fermion_indices = fermion_indices,
                                         boson_indices = boson_indices,
                                         bits_per_boson = bits_per_boson,
                                         params = {'verbosity' : 2})

# Hamiltonian at t > 0
H = H0 + g * (n('up',0) + n('dn',0)) * (a_dag('B',0) + a('B',0))

## Construct a Solver object
S = Solver(gf_struct, chi_indices, t_max = t_max, n_t = n_t)

# Set initial state
S.initial_state = init_state

# Parameters for GF/Chi calculation
gf_params = {}
gf_params['verbosity'] = 2                      # Verbosity level
gf_params['hbar'] = 1.0                         # Planck's constant
gf_params['hamiltonian_interpol'] = 'Trapezoid' # Trapezoid rule interpolation of H(t)
gf_params['lanczos_min_matrix_size'] = 33       # Use LAPACK for subspaces of dim 32 and smaller

# Run calculation!
S.compute_2t_obs(h = H, params = gf_params)

# On MPI rank 0, save results to HDF5 archive
if mpi.is_master_node():
    with HDFArchive('example.h5', 'w') as ar:
        ar['init_state'] = init_state
        ar['g_l'] = S.g_l
        ar['g_g'] = S.g_g
        ar['chi'] = S.chi
