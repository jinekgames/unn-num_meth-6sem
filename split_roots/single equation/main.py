# solution run

# import params
from unittest import TestCase
from params import EPS, FUNC, RANGE, MAX_COUNTER_VALUE, PHI
# import process func
from iteration_methods import Bisec, BasicIter
# import plot library
import matplotlib.pyplot as plot
import numpy as np

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

        x = Bisec(FUNC, RANGE, EPS, MAX_COUNTER_VALUE)

        print(
            "Root was successfully found:\n",
            "x    = ", x, " +- ", EPS, "\n",
            "f(x) = ", FUNC(x),
            sep=""
        )

    except BaseException:
        print(TextException())


    # Plot
    t = np.arange(RANGE.start, RANGE.end, 0.1)
    fig, ax = plot.subplots()
    ax.set_title("Root: " + str(x) + " +- " + str(EPS))
    plot.plot(t, FUNC(t), "r--")
    plot.plot([x], [FUNC(x)], "bo")
    plot.show()



if __name__ == "__main__":
    main()
    print("мучас грасиас офишион ешто паравасотрощ шууууууууууу")
