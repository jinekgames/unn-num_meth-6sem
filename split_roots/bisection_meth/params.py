# params for equation solution

# Imports

from math import exp

# import some structures
from base_objs import Range



# max iterations count
MAX_COUNTER_VALUE = 20

# accuracy
EPS     = 1e-4

# equation function
FUNC    = lambda x : x**2 - exp(-x)

# root search range
RANGE   = Range(0.5, 1)

# iteration function (for basic iter meth)
PHI     = lambda x : exp( -x/2 )
