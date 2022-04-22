# solution run

# import params
from unittest import TestCase
from params import *
# import process func
from solution_methods import BasicIter, Newtone
# import plot library
import matplotlib.pyplot as plot
import numpy as np

from base_objs import TextException, CalculateVectorFunction



# here the calculations start
def main():

    # check if data is correct
    if len(FUNC) != len(X0) and len(X0) != len(PHI):
        raise "Incorrect data (check params.py)"


    """
    Looking for root vector with Basic Iteration method
    """

    #try:

    x = BasicIter(FUNC, PHI, X0, EPS, MAX_COUNTER_VALUE)

    print(
        "\nRoot was successfully found:\n",
        "x    = ", x, " +- ", EPS, "\n",
        "f(x) = ", CalculateVectorFunction(FUNC, x), "\n\n",
        sep=""
    )

    #except BaseException:
    #    print(TextException())


    """
    Now Newtone method
    """

    x = Newtone(func, X0, EPS)

    print(
        "\nRoot was successfully found:\n",
        "x    = ", x, " +- ", EPS, "\n\n",
        sep=""
    )



if __name__ == "__main__":
    main()
    for i in X0:
        print("мучас грасиас офишион ешто паравасотрощ шууууууууууу")
