# methods for root refinement for f of equations

from math import fabs, sin, cos
from pickle import TRUE
import sympy as sp
import numpy as np
# setting symbols for symbol functions
x,y,z,t = sp.symbols('x y z t')

# import some structures
from ast import Lambda
from base_objs import *


# default max count of cycle runs
ITERATIONS_MAX_COUNT = 52

# Debug is ON flag
DEBUG = TRUE



def BasicIter(f: list, phi: list, x0: list, accurasy: float, iteration_max_count: int = ITERATIONS_MAX_COUNT) -> list:

    """
    Basic iteration method
    phi is iteration function
    """

    # to break unlim cycles
    counter = 0


    # accuracy validation
    assert accurasy < 1, "Accuracy should be less then 1"


    # debug
    #calculating column size
    maxIntLen = 0
    for el in x0:
        l = GetIntDigitsCount(el)
        if l > maxIntLen:
            maxIntLen = l
    vec_len = len(x0)
    tableVecColSize = (maxIntLen + GetFloatDigitsCount(accurasy) + ACCURACY_ROUND_OFFSET)*vec_len + 2*vec_len + 1
    tableScalColSize = maxIntLen + GetFloatDigitsCount(accurasy) + ACCURACY_ROUND_OFFSET + 1
    if (DEBUG):
        print("\n\nBasic Iter meth for eq. system\n")
        print("X0:\t", x0, sep="")
        print()
        print("┌" + "─"*3, "─"*tableVecColSize, "─"*tableScalColSize + "─┐", sep="─┬─")
        print(
            "│" + "k".center(3),
            "x".center(tableVecColSize),
            "Delta".center(tableScalColSize) + " │",
            sep=" │ "
        )



    x_prev = []
    x_next = x0


    # calculations themself
    while (counter < iteration_max_count):

        x_prev = x_next
        x_next = CalculateVectorFunction(phi, x_prev)

        # debug
        if (DEBUG):
            print("├" + "─"*3, "─"*tableVecColSize, "─"*tableScalColSize + "─┤", sep="─┼─")
            print(
                "│" + str(counter).rjust(3),
                str(RoundByAccurVec( x_next,                            accurasy )).center(tableVecColSize),
                str(RoundByAccur(    VectorsMaxDelta(x_next, x_prev),   accurasy )).rjust(tableScalColSize) + " │",
                sep=" │ "
            )

        # root found
        if VectorsMaxDelta(x_next, x_prev) < accurasy:
            if (DEBUG):
                print("└" + "─"*3, "─"*tableVecColSize, "─"*tableScalColSize + "─┘", sep="─┴─")
            return RoundByAccurVec(x_next, accurasy)

        counter += 1

    raise Exception(
        "Root was not found after " + str(counter) + " iterations"              + "\n" + \
        "Accuracy was not reached"                                              + "\n" + \
        "The approcsimate value: "  + str(RoundByAccurVec(x_next, accurasy))    + "\n" + \
        "Accuracy:               "  + str(VectorsMaxDelta(x_next, x_prev))
    )


def Newtone(f: list, x0: list, accurasy: float, iteration_max_count: int = ITERATIONS_MAX_COUNT) -> list:
    
    """
    Solving equation f using Newtone meth
    U should use symbol functions defined in sympy, because the function does calculate Yacoby matrix by itself
    """

    # accuracy validation
    assert accurasy < 1, "Accuracy should be less then 1"


    # debug
    #calculating column size
    maxIntLen = 0
    for el in x0:
        l = GetIntDigitsCount(el)
        if l > maxIntLen:
            maxIntLen = l
    vec_len = len(x0)
    tableVecColSize = (maxIntLen + GetFloatDigitsCount(accurasy) + ACCURACY_ROUND_OFFSET)*vec_len + 2*vec_len + 1
    tableScalColSize = maxIntLen + GetFloatDigitsCount(accurasy) + ACCURACY_ROUND_OFFSET + 1
    if (DEBUG):
        print("\n\nNewtone meth for eq. system\n")
        print("X0:\t", x0, sep="")
        print()
        print("┌" + "─"*3, "─"*tableVecColSize, "─"*tableScalColSize + "─┐", sep="─┬─")
        print(
            "│" + "k".center(3),
            "x".center(tableVecColSize),
            "Delta".center(tableScalColSize) + " │",
            sep=" │ "
        )


    jac     = [
        [ sp.diff(f[0], x),    sp.diff(f[0], y) ],
        [ sp.diff(f[1], x),    sp.diff(f[1], y) ],
    ]
    jac_num = sp.lambdify([x, y], jac)
    f       = sp.lambdify([x, y], f)
    i       = x0
    
    # to break unlim cycles
    counter = 0

    while (counter < iteration_max_count):

        diff = np.linalg.solve(
            np.array(jac_num(i[0], i[1])),
            -np.array(f(i[0], i[1]))
        )
        i += diff

        # debug
        if (DEBUG):
            print("├" + "─"*3, "─"*tableVecColSize, "─"*tableScalColSize + "─┤", sep="─┼─")
            print(
                "│" + str(counter).rjust(3),
                str(RoundByAccurVec( i,                     accurasy )).center(tableVecColSize),
                str(RoundByAccur(    np.linalg.norm(diff),  accurasy )).rjust(tableScalColSize) + " │",
                sep=" │ "
            )

        if np.linalg.norm(diff) < accurasy:
            if (DEBUG):
                print("└" + "─"*3, "─"*tableVecColSize, "─"*tableScalColSize + "─┘", sep="─┴─")
            return RoundByAccurVec(i, accurasy)

        counter += 1
