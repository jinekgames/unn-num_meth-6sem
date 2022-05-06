# Runge-Kutta methods for systems of differential equations solution

import sympy as sp
import numpy as np

# import some structures
from base_objs import *


# DEBUG lvl (0 - off,  1 - simple,  2 - full)
DEBUG = 1


# (DEBUG) default size for column with float type value
DEFAULT_FLOAT_COLUMN_SIZE          = 11
DEFAULT_FLOAT_COLUMN_SIZE_INT_PART = 3



def RungeKutta2(f: list, x0: list, h: float, n: int) -> list:

    """
    Solving differentioal equations system using 2-order Runge-Kutta method

    Get     functions vector (our system),
            Koshi task vector (default x = x0 value, y0 = y(x0), z0 = z(x0)),
            segment size
            and number of segments (number of dots - 1)

    Return  vector of point of our solution
    """


    # DEBUG
    if (DEBUG):
        print("\n\n2-order Runge-Kutta meth for dif. eq. system\n")
        print("Koshi:", x0, sep="\t")
        print("step =", h)
        print("number of steps:", n)
        print("range:  from  ", x0[0], "  to  ", x0[0] + h*n, sep="")
    if (DEBUG >= 2):
        print()
        tableFloatColCount = 3
        tableFloatColSize = DEFAULT_FLOAT_COLUMN_SIZE
        print("┌" + "─"*3, ("─┬─" + "─"*tableFloatColSize)*tableFloatColCount, "─┐", sep="")
        print(
            "│" + "k".center(3),
            "x".center(tableFloatColSize),
            "y".center(tableFloatColSize),
            "z".center(tableFloatColSize),
            sep=" │ ", end=" │\n"
        )


    # vector of point of our solution
    x = [ x0, ]


    for i in range(n):

        k1 = CalculateVectorFunction(f, x[i])
        k2 = CalculateVectorFunction(f, [ x[i][0] + h,  x[i][1] + h*k1[0],  x[i][2] + h*k1[1] ])
        
        x.append([
            x[i][0] + h,
            x[i][1] + h/2 * (k1[0] + k2[0]),
            x[i][2] + h/2 * (k1[1] + k2[1]),
        ])


        # DEBUG
        if (DEBUG >= 2):
            print("├" + "─"*3, ("─┼─" + "─"*tableFloatColSize)*tableFloatColCount, "─┤", sep="")
            print(
                "│" + str(i).rjust(3),
                str(RoundBySymbCount(x[i+1][0], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][1], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][2], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                sep=" │ ", end=" │\n"
            )


    if (DEBUG >= 2):
        print("└" + "─"*3, ("─┴─" + "─"*tableFloatColSize)*tableFloatColCount, "─┘", sep="")


    return x


def RungeKutta3(f: list, x0: list, h: float, n: int) -> list:

    """
    Solving differentioal equations system using 3-order Runge-Kutta method

    Get     functions vector (our system),
            Koshi task vector (default x = x0 value, y0 = y(x0), z0 = z(x0)),
            segment size
            and number of segments (number of dots - 1)

    Return  vector of point of our solution
    """


    # DEBUG
    if (DEBUG):
        print("\n\n3-order Runge-Kutta meth for dif. eq. system\n")
        print("Koshi:", x0, sep="\t")
        print("step =", h)
        print("number of steps:", n)
        print("range:  from  ", x0[0], "  to  ", x0[0] + h*n, sep="")
    if (DEBUG >= 2):
        print()
        tableFloatColCount = 3
        tableFloatColSize = DEFAULT_FLOAT_COLUMN_SIZE
        print("┌" + "─"*3, ("─┬─" + "─"*tableFloatColSize)*tableFloatColCount, "─┐", sep="")
        print(
            "│" + "k".center(3),
            "x".center(tableFloatColSize),
            "y".center(tableFloatColSize),
            "z".center(tableFloatColSize),
            sep=" │ ", end=" │\n"
        )


    # vector of point of our solution
    x = [ x0, ]


    for i in range(n):

        k1 = CalculateVectorFunction(f, x[i])
        k2 = CalculateVectorFunction(f, [ x[i][0] + h/2,  x[i][1] + h/2*k1[0],            x[i][2] + h/2*k1[1] ])
        k3 = CalculateVectorFunction(f, [ x[i][0] + h,    x[i][1] - h*k1[0] + 2*h*k2[0],  x[i][2] - h*k1[1] + 2*h*k2[1] ])
        
        x.append([
            x[i][0] + h,
            x[i][1] + h/6 * (k1[0] + 4*k2[0] + k3[0]),
            x[i][2] + h/6 * (k1[1] + 4*k2[1] + k3[1]),
        ])


        # DEBUG
        if (DEBUG >= 2):
            print("├" + "─"*3, ("─┼─" + "─"*tableFloatColSize)*tableFloatColCount, "─┤", sep="")
            print(
                "│" + str(i).rjust(3),
                str(RoundBySymbCount(x[i+1][0], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][1], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][2], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                sep=" │ ", end=" │\n"
            )


    if (DEBUG >= 2):
        print("└" + "─"*3, ("─┴─" + "─"*tableFloatColSize)*tableFloatColCount, "─┘", sep="")


    return x


def RungeKutta4(f: list, x0: list, h: float, n: int) -> list:

    """
    Solving differentioal equations system using 4-order Runge-Kutta method

    Get     functions vector (our system),
            Koshi task vector (default x = x0 value, y0 = y(x0), z0 = z(x0)),
            segment size
            and number of segments (number of dots - 1)

    Return  vector of point of our solution
    """


    # DEBUG
    if (DEBUG):
        print("\n\n4-order Runge-Kutta meth for dif. eq. system\n")
        print("Koshi:", x0, sep="\t")
        print("step =", h)
        print("number of steps:", n)
        print("range:  from  ", x0[0], "  to  ", x0[0] + h*n, sep="")
    if (DEBUG >= 2):
        print()
        tableFloatColCount = 3
        tableFloatColSize = DEFAULT_FLOAT_COLUMN_SIZE
        print("┌" + "─"*3, ("─┬─" + "─"*tableFloatColSize)*tableFloatColCount, "─┐", sep="")
        print(
            "│" + "k".center(3),
            "x".center(tableFloatColSize),
            "y".center(tableFloatColSize),
            "z".center(tableFloatColSize),
            sep=" │ ", end=" │\n"
        )


    # vector of point of our solution
    x = [ x0, ]


    for i in range(n):

        k1 = CalculateVectorFunction(f, x[i])
        k2 = CalculateVectorFunction(f, [ x[i][0] + h/2,  x[i][1] + h/2*k1[0],  x[i][2] + h/2*k1[1] ])
        k3 = CalculateVectorFunction(f, [ x[i][0] + h/2,  x[i][1] + h/2*k2[0],  x[i][2] + h/2*k2[1] ])
        k4 = CalculateVectorFunction(f, [ x[i][0] + h,    x[i][1] + h*k3[0],    x[i][2] + h*k3[1] ])
        
        x.append([
            x[i][0] + h,
            x[i][1] + h/6 * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]),
            x[i][2] + h/6 * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]),
        ])


        # DEBUG
        if (DEBUG >= 2):
            print("├" + "─"*3, ("─┼─" + "─"*tableFloatColSize)*tableFloatColCount, "─┤", sep="")
            print(
                "│" + str(i).rjust(3),
                str(RoundBySymbCount(x[i+1][0], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][1], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][2], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                sep=" │ ", end=" │\n"
            )


    if (DEBUG >= 2):
        print("└" + "─"*3, ("─┴─" + "─"*tableFloatColSize)*tableFloatColCount, "─┘", sep="")


    return x

