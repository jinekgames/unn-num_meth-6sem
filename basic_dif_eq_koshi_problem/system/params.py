"""
Params for sample solution
"""


import numpy as np
import sympy as sp
from base_objs import *



# equation functions vector
FUNC    = [
    lambda x :   x[2] * x[2] + x[0] ,       # y' = z^2 + x;
    lambda x :   x[0] * x[1]        ,       # z' = xy
]

# search range
RANGE = Range(0, 1)

# Koshi condition
X0 = [ 0, 1, 0.5 ]      # x0, y0, z0

SEGMENTS_COUNT = 10

