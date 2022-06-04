# params for equation solution

# Imports

from math import exp
import numpy as np
import sympy as sp

# import some structures
from base_objs import Range



# accuracy
EPS     = 1e-4

# equation functions vector
FUNC    = [
    lambda x :   np.cos(x[1])     + 2 * x[0]   - 2   ,
    lambda x :   np.cos(x[0] - 1) + x[1]       - 1.6 ,
]

from solution_methods import Newtone, x, y
func = [
    sp.cos(y)     + 2 * x   - 2   ,
    sp.cos(x - 1) + y       - 1.6 ,
]

"""
resolution:
x = -0.1335583261035361
y = -0.8665808075256101
"""

# starting search point
X0      = [ 100, 100, ]

# iteration functions vector (for basic iter meth)
PHI     = [
    lambda x : (-np.cos(x[1]) + 2) / 2   ,
    lambda x : -np.cos(x[0] - 1) + 1.6  ,
]
