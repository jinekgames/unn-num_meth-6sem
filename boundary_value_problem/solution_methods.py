"""
here r some function whick solve some math tasks
"""


import numpy as np

from base_objs import *



DEBUG = FALSE

MAX_ITERATIONS_COUNT = 52



def LinAlgRelaxIter(x: list, a: list, b: list, w: float, n: float) -> list:

    """
    Making iteration of solving lin alg system using relaxation method

    Arguements:

         x:  array  (n)
        previous step of solution

         a:  array of array  (n*n)
        matrix of system
        should be square

         b:  array  (n)
        right column

         w:  float
        some const
        w = (0, ..., 1) -- Lower Relaxation
        w = 1           -- Zeydel method
        max_w = 2

         n:  float
        size of arrays

    Return: 

         array  (n)
        next step of solution

    """


    sol = np.zeros(n)

    # iteration formula itself
    for i in range(n):

        s1 = 0
        for j in range(i):
            s1 += a[j][i] * sol[j]

        s2 = 0
        for j in range(i + 1, n):
            s2 += a[j][i] * x[j]

        sol[i] = ( -w * s1  +  (1 - w) * a[i][i] * x[i]  -  w * s2  +  w * b[i] ) / a[i][i]


    return sol


def LinAlgRelax(x: list, a: list, b: list = 0, accuracy: float = 1e-5, w: float = 1) -> list:

    """
    Solving lin alg system using relaxation iteration method

    Arguements:

         x:  array  (n)
        previous step of solution

         a:  array of array  (n*n)
        matrix of system
        should be square

         b:  array  (n)
        right column

         accuracy: float
        accuracy of calculations

         w:  float
        some const
        w = (0, ..., 1) -- Lower Relaxation
        w = 1           -- Zeydel method
        max_w = 2

    Return: 

         array  (n)
        next step of solution

    """


    n = len(x)

    if (not len(b)):
        b = np.zeros(n)

    assert (n == len(a) and n == len(a[0]) and n == len(b)),  "incorrect data sizes"
    assert (w > 0 and w < 2),  "incorrect w"

    x_next = np.zeros(n)
    x_prev = x

    counter = 0
    while TRUE:   # many iterations for potryasayushaya accuracy
        x_next = LinAlgRelaxIter(x_prev, a, b, w, n)
        if (VectorsMaxDelta(x_next, x_prev) < accuracy):
            return RoundByAccurVec(x_next, accuracy)
        if (counter >= MAX_ITERATIONS_COUNT):
            return RoundByAccurVec(x_next, x_next - x_prev)
        x_prev = x_next


def BoundValPuassonRectEq(f: lambda_type, rect: RECT_DOUBLE, conds: RECT_FUNC, x_step: float, y_step = 0) -> list:

    """
    This function can calculate Puasson equation using numerical grid method
    Linear alg system is being calculated using LowRelaxation methon (TODO: not optimized)

    Arguements:

         f:       lambda
        function in right side of Puasson eq.                left  right  
                                                               +——————— >  X
         rect:    RECT of floats                        bottom | ┌─┐
        this value shows edges of area                     top | └─┘
                                                            Y  V        // carthesian system is computer like (OY goes down)             
         conds:   RECT of lambdas
        conditions on the rect edges

         x_step:  float
        size of grid on Ox

         y_step:  float
        size of grid on Oy
        if was not set then = x_step

        x_step and y_step must be such that we'll have the same count of segments on Oy and Ox axis

    Return:

         array
        the array with counted values on points of our grid inside the rect
        they are numbered from 0 to i_max and does upside down according to the picture, from left to right

    """         


    # check edges and steps
    assert int((rect.right - rect.left) / x_step) == ((rect.right - rect.left) / x_step),  "incorrect x_step"
    assert int((rect.top - rect.bottom) / y_step) == ((rect.top - rect.bottom) / y_step),  "incorrect y_step"

    # calculate max values for j, k and i
    n = int( (rect.right - rect.left) / x_step ) + 1
    m = int( (rect.top - rect.bottom) / y_step ) + 1
    i_max = (n - 2) * (m - 2)
    print("n =", n, "m =", m, "i_max =", i_max, end="\n\n")
    assert (n == m),  "incorrect grid"

    # macros to calculate one coords from another
    j2x  =  lambda j :     j * x_step
    k2y  =  lambda k :     k * y_step
    kj2i =  lambda k, j : (k - 1) * (n - 2) + j
    i2jk =  lambda i:    ( i % (n - 2) + 1,  int(i / (n - 2)) + 1 )

    # define matrix of lin eq sys
    matr = np.zeros((i_max, i_max))
    rvec = np.zeros(i_max)
    

    # fill the matrix
    for i in range(i_max):

        j, k = i2jk(i)


        # f_i
        rvec[i] += f([j2x(j), k2y(k)])

        # u_j_k
        matr[i][i] = - 2 / x_step**2 - 2 / y_step**2

        # u_j+1_k
        if (j + 1 == n - 1):
            rvec[i] -= conds.right([j2x(j), k2y(k)]) * 1/x_step**2
        else:
            matr[i][i+1] = 1 / x_step**2

        # u_j-1_k
        if (j - 1 == 0):
            rvec[i] -= conds.bottom([j2x(j), k2y(k)]) * 1/x_step**2
        else:
            matr[i][i-1] = 1 / x_step**2

        # u_j_k+1
        if (k + 1 == m - 1):
            rvec[i] -= conds.top([j2x(j), k2y(k)]) * 1/y_step**2
        else:
            matr[i][i + (n-2)] = 1 / y_step**2

        # u_j_k-1
        if (k - 1 == 0):
            rvec[i] -= conds.left([j2x(j), k2y(k)]) * 1/y_step**2
        else:
            matr[i][i - (n-2)] = 1 / y_step**2


    # print the matrix
    if (DEBUG):
        print("\nThe matrix of linear equation system:\n")
        for k in range(i_max):
            for j in range(i_max):
                if(matr[k][j]):
                    print(str(round(matr[k][j], 2)).rjust(6), end="")
                else:
                    print(" "*6, end="")
            print(" |", str(round(rvec[k], 2)).rjust(6), sep="")


    """
    Solve the system
    """
    x = np.array([1 for i in range(i_max)])
    x = LinAlgRelax(x, matr, rvec, 1e-6, 0.9)

    return x
