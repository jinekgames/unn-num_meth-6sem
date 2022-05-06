# Runge-Kutta methods for systems of differential equations

from types import LambdaType
import sympy as sp
import numpy as np

# import some structures
from base_objs import *


# DEBUG lvl (0 - off,  1 - simple,  2 - full)
DEBUG = 2


# (DEBUG) default size for column with float type value
DEFAULT_FLOAT_COLUMN_SIZE          = 11
DEFAULT_FLOAT_COLUMN_SIZE_INT_PART = 3



def RungeKutta2(f: lambda_type, x0: list, h: float, n: int) -> list:

    """
    Solving differentioal equations using 2-order Runge-Kutta method

    Get     function ( y' = f(list : [x, y]) ),
            Koshi task vector (default x = x0 value, y0 = y(x0)),
            segment size
            and number of segments (number of dots - 1)

    Return  vector of points of our solution
    """


    # DEBUG
    if (DEBUG):
        print("\n\n2-order Runge-Kutta meth for dif. eq.\n")
        print("Koshi:", x0, sep="\t")
        print("step =", h)
        print("number of steps:", n)
        print("range:  from  ", x0[0], "  to  ", x0[0] + h*n, sep="")
    if (DEBUG >= 2):
        print()
        tableFloatColCount = 2
        tableFloatColSize = DEFAULT_FLOAT_COLUMN_SIZE
        print("┌" + "─"*3, ("─┬─" + "─"*tableFloatColSize)*tableFloatColCount, "─┐", sep="")
        print(
            "│" + "k".center(3),
            "x".center(tableFloatColSize),
            "y".center(tableFloatColSize),
            sep=" │ ", end=" │\n"
        )


    # vector of points of our solution
    x = [ x0, ]


    for i in range(n):

        k1 = f(x[i])
        k2 = f([ x[i][0] + h,  x[i][1] + h*f(x[i]) ])

        
        x.append([
            x[i][0] + h,
            x[i][1] + h/2 * (k1 + k2)
        ])


        # DEBUG
        if (DEBUG >= 2):
            print("├" + "─"*3, ("─┼─" + "─"*tableFloatColSize)*tableFloatColCount, "─┤", sep="")
            print(
                "│" + str(i).rjust(3),
                str(RoundBySymbCount(x[i+1][0], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][1], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                sep=" │ ", end=" │\n"
            )


    if (DEBUG >= 2):
        print("└" + "─"*3, ("─┴─" + "─"*tableFloatColSize)*tableFloatColCount, "─┘", sep="")


    return x


def RungeKutta3(f: lambda_type, x0: list, h: float, n: int) -> list:

    """
    Solving differentioal equations using 3-order Runge-Kutta method

    Get     function ( y' = f(list : [x, y]) ),
            Koshi task vector (default x = x0 value, y0 = y(x0)),
            segment size
            and number of segments (number of dots - 1)

    Return  vector of points of our solution
    """


    # DEBUG
    if (DEBUG):
        print("\n\n3-order Runge-Kutta meth for dif. eq.\n")
        print("Koshi:", x0, sep="\t")
        print("step =", h)
        print("number of steps:", n)
        print("range:  from  ", x0[0], "  to  ", x0[0] + h*n, sep="")
    if (DEBUG >= 2):
        print()
        tableFloatColCount = 2
        tableFloatColSize = DEFAULT_FLOAT_COLUMN_SIZE
        print("┌" + "─"*3, ("─┬─" + "─"*tableFloatColSize)*tableFloatColCount, "─┐", sep="")
        print(
            "│" + "k".center(3),
            "x".center(tableFloatColSize),
            "y".center(tableFloatColSize),
            sep=" │ ", end=" │\n"
        )


    # vector of points of our solution
    x = [ x0, ]


    for i in range(n):

        k1 = f(x[i])
        k2 = f([ x[i][0] + h/2,  x[i][1] + h/2*k1 ])
        k3 = f([ x[i][0] + h,    x[i][1] - h*k1 + 2*h*k2 ])

        
        x.append([
            x[i][0] + h,
            x[i][1] + h/6 * (k1 + 4*k2 + k3)
        ])


        # DEBUG
        if (DEBUG >= 2):
            print("├" + "─"*3, ("─┼─" + "─"*tableFloatColSize)*tableFloatColCount, "─┤", sep="")
            print(
                "│" + str(i).rjust(3),
                str(RoundBySymbCount(x[i+1][0], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][1], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                sep=" │ ", end=" │\n"
            )


    if (DEBUG >= 2):
        print("└" + "─"*3, ("─┴─" + "─"*tableFloatColSize)*tableFloatColCount, "─┘", sep="")


    return x


def RungeKutta4(f: lambda_type, x0: list, h: float, n: int) -> list:

    """
    Solving differentioal equations using 4-order Runge-Kutta method

    Get     function ( y' = f(list : [x, y]) ),
            Koshi task vector (default x = x0 value, y0 = y(x0)),
            segment size
            and number of segments (number of dots - 1)

    Return  vector of points of our solution
    """


    # DEBUG
    if (DEBUG):
        print("\n\n4-order Runge-Kutta meth for dif. eq.\n")
        print("Koshi:", x0, sep="\t")
        print("step =", h)
        print("number of steps:", n)
        print("range:  from  ", x0[0], "  to  ", x0[0] + h*n, sep="")
    if (DEBUG >= 2):
        print()
        tableFloatColCount = 2
        tableFloatColSize = DEFAULT_FLOAT_COLUMN_SIZE
        print("┌" + "─"*3, ("─┬─" + "─"*tableFloatColSize)*tableFloatColCount, "─┐", sep="")
        print(
            "│" + "k".center(3),
            "x".center(tableFloatColSize),
            "y".center(tableFloatColSize),
            sep=" │ ", end=" │\n"
        )


    # vector of points of our solution
    x = [ x0, ]


    for i in range(n):

        k1 = f(x[i])
        k2 = f([ x[i][0] + h/2,  x[i][1] + h/2*k1 ])
        k3 = f([ x[i][0] + h/2,  x[i][1] + h/2*k2 ])
        k4 = f([ x[i][0] + h,    x[i][1] + h*k3 ])

        
        x.append([
            x[i][0] + h,
            x[i][1] + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        ])


        # DEBUG
        if (DEBUG >= 2):
            print("├" + "─"*3, ("─┼─" + "─"*tableFloatColSize)*tableFloatColCount, "─┤", sep="")
            print(
                "│" + str(i).rjust(3),
                str(RoundBySymbCount(x[i+1][0], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                str(RoundBySymbCount(x[i+1][1], DEFAULT_FLOAT_COLUMN_SIZE - DEFAULT_FLOAT_COLUMN_SIZE_INT_PART - 1)).center(tableFloatColSize),
                sep=" │ ", end=" │\n"
            )


    if (DEBUG >= 2):
        print("└" + "─"*3, ("─┴─" + "─"*tableFloatColSize)*tableFloatColCount, "─┘", sep="")


    return x
