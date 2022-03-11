# iteration methods for root refinement

from math import fabs

# import some structures
from ast import Lambda
from base_objs import Range


# default max count of cycle runs
ITERATIONS_MAX_COUNT = 10

# Debug is ON flag
DEBUG = False


def GetFloatDigitsCount(x: float) -> float:
    """
    Returns count of digits after point
    """
    return len(str(x).split('.')[1])

def RoundByAccur(x: float, accuracy: float) -> float:
    """
    Rounds the float X according to given accuracy
    """
    return float(f"{x:.{GetFloatDigitsCount(accuracy)}f}")



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
    assert(a >= b, "Incorrect range: a < b")


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
            raise Exception("Root is on the edge, this case have not been designed")
        # no root in the range
        else:
            raise Exception("No roots are there in the " + range.to_str())


        # root found
        if b - a <= accurasy * 2:
            return RoundByAccur(c, accurasy)

        counter += 1

    raise Exception(
        "Root was not found after " + str(counter) + " iterations"      + "\n" + \
        "Accuracy was not reached"                                      + "\n" + \
        "The approcsimate value: "  + str(RoundByAccur(c, accurasy))    + "\n" + \
        "Accuracy:               "  + str(fabs(b - a))
    )

def BasicIter(f: Lambda, phi: Lambda, range: Range, accurasy: float, iteration_max_count: int = ITERATIONS_MAX_COUNT) -> float:

    """
    Basic iteration method
    phi is iteration function
    """

    # to break unlim cycles
    counter = 0

    # range validation
    assert(range.start >= range.end, "Incorrect range: a < b")


    # debug
    if (DEBUG):
        print("Start of range\t:\t", range.start, sep="")
        print("End of range\t:\t", range.end, sep="")


    x_prev = 0
    x_next = range.start


    # calculations themself
    while (counter < iteration_max_count):

        x_prev = x_next
        x_next = phi(x_prev)

        # debug
        if (DEBUG):
            print("x_k =",   x_prev)
            print("x_k+1 =", x_next)

        # out of range
        assert(x_next <= range.end, "No roots are there in the " + range.to_str())

        # root found
        if fabs(x_next - x_prev) < accurasy:
            return RoundByAccur(x_next, accurasy)

        counter += 1

    raise Exception(
        "Root was not found after " + str(counter) + " iterations"          + "\n" + \
        "Accuracy was not reached"                                          + "\n" + \
        "The approcsimate value: "  + str(RoundByAccur(x_next, accurasy))    + "\n" + \
        "Accuracy:               "  + str(fabs(x_next - x_prev))
    )
