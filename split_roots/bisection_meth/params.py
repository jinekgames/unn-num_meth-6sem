# params for equation solution

# Imports

from math import exp

# import some structures
from base_objs import Range



# max iterations count
MAX_COUNTER_VALUE = 20

# accuracy
EPS     = 1e-2

# equation function
FUNC    = lambda x : x**5 + 2*x - 8

# root search range
RANGE   = Range(0, 3)

# iteration function (for basic iter meth)
PHI     = lambda x : exp( -x/2 )
