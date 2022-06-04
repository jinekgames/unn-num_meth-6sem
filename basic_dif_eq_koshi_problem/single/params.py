"""
Params for sample solution
"""


import numpy as np
import sympy as sp
from base_objs import *



# equation functions vector
FUNC = lambda x :   2 * x[1] / (x[0] * np.log(x[0]))  +  1 / x[0]
# solution range
RANGE = Range(2, 3)

# Koshi condition
X0 = [ 2, -np.log(2) ]      # x0, y0 = y(x0)

SEGMENTS_COUNT = 10
