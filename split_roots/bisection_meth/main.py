# solution run

# import params
from unittest import TestCase
from params import EPS, FUNC, RANGE, MAX_COUNTER_VALUE, PHI
# import process func
from iteration_methods import Bisec, BasicIter
# import plot library
#from matplotlib import pyplot

import sys
import linecache



# exception text compiler
def TextException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    return 'EXCEPTION IN ({}, LINE {} "{}"):\n{}'.format(filename, lineno, line.strip(), exc_obj)



# here the calculations start
def main():

    try:

        x = BasicIter(FUNC, PHI, RANGE, EPS, MAX_COUNTER_VALUE)

        print(
            "Root was successfully found:\n",
            "x = ", x, "+-", EPS,
            sep=""
        )

    except BaseException:
        print(TextException())


    # TODO: plot



if __name__ == "__main__":
    main()
