# params for equation solution

# Imports

from math import exp, log10, log, fabs
import numpy as np

# import some structures
from base_objs import Range



# accuracy
EPS      = 1e-1
EPS_SPEC = 1e-3

# root search range
RANGE   = Range(0.1, 10)

# equation function
FUNC       = lambda x : x + np.log10(x) - 0.5
FUNC_DERIV = lambda x : 1 + 1 / (log(10)*x)

# iteration function (for basic iter meth)
PHI        = lambda x : 0.5 - np.log10(x)
