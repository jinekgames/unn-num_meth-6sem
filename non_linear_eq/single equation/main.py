# solution run

# import params
from unittest import TestCase

from sympy import false
from params import *
# import process func
from iteration_methods import *
# import plot library
import matplotlib.pyplot as plot
import numpy as np

#from ..common.objects import TextException
from base_objs import TextException


# set true if u wanna see plots
SHOW_PLOTS = false



# here the calculations start
def main():

    """
    At first we are looking for approximate value of the root using bisection method
    """

    try:

        x = Bisec(FUNC, RANGE, EPS)

        print(
            "\nRoot was successfully found:\n",
            "x    = ", x, " +- ", EPS, "\n",
            "f(x) = ", FUNC(x), "\n\n",
            sep=""
        )

    except BaseException:
        print(TextException())

    # Plot
    if (SHOW_PLOTS):
        ShowRootPlot(FUNC, x, RANGE, EPS)

    

    """
    Now we will specify our root using basic iteration method and the value we got in the first part
    """

    try:

        x_snd = BasicIter(FUNC, PHI, x, EPS_SPEC)

        print(
            "\nRoot was successfully found:\n",
            "x    = ", x_snd, " +- ", EPS_SPEC, "\n",
            "f(x) = ", FUNC(x_snd), "\n\n",
            sep=""
        )

    except BaseException:
        print(TextException())

    # Plot
    if (SHOW_PLOTS):
        ShowRootPlot(FUNC, x_snd, RANGE, EPS_SPEC)

    

    """
    And now we will also specify the root. But now we'll use Newtone method
    """

    try:

        x_trd = Newtone(FUNC, FUNC_DERIV, x, EPS_SPEC)

        print(
            "\nRoot was successfully found:\n",
            "x    = ", x_trd, " +- ", EPS_SPEC, "\n",
            "f(x) = ", FUNC(x_snd), "\n\n",
            sep=""
        )

    except BaseException:
        print(TextException())

    # Plot
    if (SHOW_PLOTS):
        ShowRootPlot(FUNC, x_snd, RANGE, EPS_SPEC)




if __name__ == "__main__":
    main()
    print("\n"*10 + "мучас грасиас офишион ешто паравасотрощ шууууууууууу")
