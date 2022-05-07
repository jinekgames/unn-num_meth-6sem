# iteration methods for root refinement for single equation

from math import fabs

# import some structures
from ast import Lambda
from pickle import TRUE

from sympy import false, true
#from ..common.objects import *
from base_objs import *


# default max count of cycle runs
ITERATIONS_MAX_COUNT = 52
# float type contains 52 bytes of lyterals
# thats why we can devide 52 times to get the maximum accuracy in bython float type

# Debug is ON flag
DEBUG = TRUE

# Count of space symbols in the table Debug output
TABLE_INT_PART_SIZE = 5



def Bisec(f: Lambda, range: Range, accurasy: float, iteration_max_count: int = ITERATIONS_MAX_COUNT) -> float:

    """
    Bisection method
    """

    # to break unlim cycles
    counter = 0

    # searching range
    a = range.start
    b = range.end

    # range validation
    assert a <= b, "Incorrect range: a < b"

    # accuracy validation
    assert accurasy < 1, "Accuracy should be less then 1"


    # debug
    tableColSize = GetFloatDigitsCount(accurasy) + GetIntDigitsCount(f((range.start + range.end)/2)) + 1 + ACCURACY_ROUND_OFFSET
    if (DEBUG):
        print("\n\nBisection meth for eq.\n")
        print("Start of range\t:\t", range.start, sep="")
        print("End of range\t:\t", range.end, sep="")
        print()
        print("┌" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┐", sep="─┬─")
        print(
            "│" + "k".center(3),
            "x".center(tableColSize),
            "f(x)".center(tableColSize),
            "Delta".center(tableColSize) + " │",
            sep=" │ "
        )


    # calculations themself
    while (counter < iteration_max_count):

        c = (a + b) / 2

        # root is locatied in the left part
        if   f(a) * f(c) < 0:
            b = c
        # root is located in the right part
        elif f(b) * f(c) < 0:
            a = c
        # one of the range edges if the root
        elif f(a) == 0 or f(b) == 0:
            raise Exception("Root is on the edge, this case have not been designed yet")
        # no root in the range
        else:
            raise Exception("No roots are there in the " + range.to_str())

            
        # debug
        if (DEBUG):
            print("├" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┤", sep="─┼─")
            print(
                "│" + str(counter).rjust(3),
                str(RoundByAccur( c,                 accurasy )).rjust(tableColSize),
                str(RoundByAccur( f(c),              accurasy )).rjust(tableColSize),
                str(RoundByAccur( fabs((b) - (a)),   accurasy )).rjust(tableColSize) + " │",
                sep=" │ "
            )


        # root found
        if fabs((b) - (a)) <= accurasy:
            # debug
            if (DEBUG):
                print("└" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┘", sep="─┴─")

            return RoundByAccur(c, accurasy)

        counter += 1

    raise Exception(
        "Root was not found after " + str(counter) + " iterations"      + "\n" + \
        "Accuracy was not reached"                                      + "\n" + \
        "The approcsimate value: "  + str(RoundByAccur(c, accurasy))    + "\n" + \
        "Accuracy:               "  + str(fabs(b - a))
    )

def BasicIter(f: Lambda, phi: Lambda, x_start: float, accurasy: float, iteration_max_count: int = ITERATIONS_MAX_COUNT) -> float:

    """
    Basic iteration method
    phi is iteration function
    """

    # to break unlim cycles
    counter = 0

    # accuracy validation
    assert accurasy < 1, "Accuracy should be less then 1"


    # debug
    tableColSize = GetFloatDigitsCount(accurasy) + GetIntDigitsCount(f(x_start)) + 1 + ACCURACY_ROUND_OFFSET
    if (DEBUG):
        print("\n\nBasic Iter meth for eq.\n")
        print("Start of search\t:\t", x_start, sep="")
        print()
        print("┌" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┐", sep="─┬─")
        print(
            "│" + "k".center(3),
            "x".center(tableColSize),
            "f(x)".center(tableColSize),
            "Delta".center(tableColSize) + " │",
            sep=" │ "
        )


    x_prev = 0
    x_next = x_start


    # calculations themself
    while (counter < iteration_max_count):

        x_prev = x_next
        x_next = phi(x_prev)

        # debug
        if (DEBUG):
            print("├" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┤", sep="─┼─")
            print(
                "│" + str(counter).rjust(3),
                str(RoundByAccur( x_next,                       accurasy )).rjust(tableColSize),
                str(RoundByAccur( f(x_next),                    accurasy )).rjust(tableColSize),
                str(RoundByAccur( fabs((x_next) - (x_prev)),    accurasy )).rjust(tableColSize) + " │",
                sep=" │ "
            )
            
        # root found
        if fabs((x_next) - (x_prev)) < accurasy:
            # debug
            if (DEBUG):
                print("└" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┘", sep="─┴─")

            return RoundByAccur(x_next, accurasy)

        counter += 1

    raise Exception(
        "Root was not found after " + str(counter) + " iterations"          + "\n" + \
        "Accuracy was not reached"                                          + "\n" + \
        "The approcsimate value: "  + str(RoundByAccur(x_next, accurasy))   + "\n" + \
        "Accuracy:               "  + str(fabs(x_next - x_prev))
    )

def Newtone(f: Lambda, deriv: Lambda, x_start: float, accurasy: float, iteration_max_count: int = ITERATIONS_MAX_COUNT) -> float:

    """
    Neqtone method
    deriv is deriviative of f function
    """

    # to break unlim cycles
    counter = 0

    # accuracy validation
    assert accurasy < 1, "Accuracy should be less then 1"


    # debug
    tableColSize = GetFloatDigitsCount(accurasy) + GetIntDigitsCount(f(x_start)) + 1 + ACCURACY_ROUND_OFFSET
    if (DEBUG):
        print("\n\nNewtone meth for eq.\n")
        print("Start of search\t:\t", x_start, sep="")
        print()
        print("┌" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┐", sep="─┬─")
        print(
            "│" + "k".center(3),
            "x".center(tableColSize),
            "f(x)".center(tableColSize),
            "Delta".center(tableColSize) + " │",
            sep=" │ "
        )


    x_prev = x_start
    x_next = 0


    # calculations themself
    while (counter < iteration_max_count):

        x_next = x_prev - f(x_prev) / deriv(x_prev)

        # debug
        if (DEBUG):
            print("├" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┤", sep="─┼─")
            print(
                "│" + str(counter).rjust(3),
                str(RoundByAccur(      x_next,             accurasy )).rjust(tableColSize),
                str(RoundByAccur(    f(x_next),            accurasy )).rjust(tableColSize),
                str(RoundByAccur( fabs(x_next) - (x_prev), accurasy )).rjust(tableColSize) + " │",
                sep=" │ "
            )

        # root found
        if fabs(x_next - x_prev) < accurasy:
            # debug
            if (DEBUG):
                print("└" + "─"*3, "─"*tableColSize, "─"*tableColSize, "─"*tableColSize + "─┘", sep="─┴─")

            return RoundByAccur(x_next, accurasy)

        x_prev = x_next
        counter += 1

    raise Exception(
        "Root was not found after " + str(counter) + " iterations"          + "\n" + \
        "Accuracy was not reached"                                          + "\n" + \
        "The approcsimate value: "  + str(RoundByAccur(x_next, accurasy))   + "\n" + \
        "Accuracy:               "  + str(fabs(x_next - x_prev))
    )\
