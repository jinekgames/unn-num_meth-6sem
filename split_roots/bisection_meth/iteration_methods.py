# iteration methods for root refinement

from math import fabs

# import some structures
from ast import Lambda
from base_objs import Range


# default max count of cycle runs
ITERATIONS_MAX_COUNT = 10



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
            return c

        counter += 1

    raise Exception("Root was not found after " + counter + " iterations")


def BasicIter(f: Lambda, phi: Lambda, range: Range, accurasy: float, iteration_max_count: int = ITERATIONS_MAX_COUNT) -> float:

    """
    Basic iteration method
    phi is iteration function
    """

    # to break unlim cycles
    counter = 0

    # range validation
    assert(range.start >= range.end, "Incorrect range: a < b")


    x_prev = range.start
    x_next = 0

    # calculations themself
    while (counter < iteration_max_count):

        x_next = phi(x_prev)

        # out of range
        assert(x_next <= range.end, "No roots are there in the " + range.to_str())

        # root found
        if fabs(x_next - x_prev):
            return x_next

        counter += 1

    raise Exception("Root was not found after " + counter + " iterations")
