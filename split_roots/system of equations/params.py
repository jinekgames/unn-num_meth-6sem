# params for equation solution

# Imports

from math import exp
import numpy as np
import sympy as sp

# import some structures
from base_objs import Range



# max iterations count
MAX_COUNTER_VALUE = 20

# accuracy
EPS     = 1e-4

# equation functions vector
FUNC    = [
    lambda x :   np.cos(x[1] + 0.5) + x[0]   - 0.8 ,
    lambda x :   np.sin(x[0])       - 2*x[1] - 1.6 ,
]

from solution_methods import Newtone, x, y
func = [
    sp.cos(y + 0.5) + x   - 0.8 ,
    sp.sin(x)       - 2*y - 1.6 ,
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
    lambda x : 0.8 - np.cos(x[1] + 0.5)   ,
    lambda x : (-1.6 + np.sin(x[0])) / 2  ,
]
