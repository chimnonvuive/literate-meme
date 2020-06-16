
from matplotlib import rc

from matplotlib.pyplot import *
from numpy import *
from control import *
from scipy.signal import residue, tf2ss
from scipy.linalg import *
from sympy.solvers.solveset import linsolve
from sympy.physics.mechanics import dynamicsymbols, init_vprinting
from sympy import Matrix, symbols, sin, cos, tan, asin, acos, atan, simplify
from scipy.integrate import solve_ivp

rc('text', usetex=True)
init_vprinting()
