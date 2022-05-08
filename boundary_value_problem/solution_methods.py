"""
here r some function whick solve some math tasks
"""


import numpy as np
import scipy.linalg as lalg
from sympy import solve

from base_objs import *



DEBUG = TRUE

MAX_ITERATIONS_COUNT = 1e3



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

    if (DEBUG):
        # get exact solutuon using SciPy
        sp_x = lalg.solve(a, b.reshape((n, 1)))
        print("LinAlg solution")

    if (not len(b)):
        b = np.zeros(n)

    assert (n == len(a) and n == len(a[0]) and n == len(b)),  "incorrect data sizes"
    assert (w > 0 and w < 2),  "incorrect w"

    x_next = np.zeros(n)
    x_prev = x

    counter = 0
    while TRUE:   # many iterations for potryasayushaya accuracy
        # print(x_prev)    # ULTA SUPER MEGA DEBUG
        x_next = LinAlgRelaxIter(x_prev, a, b, w, n)
        if (VectorsMaxDelta(x_next, x_prev) < accuracy):
            if (DEBUG):
                print("Lin eq sys accuracy:", accuracy)
                print(counter + 1, "iterations done")
                print("Delta to scipy solution:", VectorsMaxDelta(sp_x, x_next))
            return RoundByAccurVec(x_next, accuracy)
        if (counter >= MAX_ITERATIONS_COUNT):
            if (DEBUG):
                print("Lin eq sys accuracy (user accuracy was not archived):", VectorsMaxDelta(sp_x, x_next))
                print(counter + 1, "iterations done")
                print("Delta to scipy solution:", VectorsMaxDelta(sp_x, x_next))
            return RoundByAccurVec(x_next, x_next - x_prev)
        x_prev = x_next
        counter += 1


def BoundValPuassonRectEq(f: lambda_type, rect: RECT_DOUBLE, conds: RECT_FUNC, x_step: float, y_step = 0, accuracy = 1e-5) -> list:

    """
    This function can calculate Puasson equation using numerical grid method
    Linear alg system is being calculated using LowRelaxation methon (TODO: add optimal w calculation)

    Arguements:

         f:        lambda
        function in right side of Puasson eq.            left  right  
                                                           +——————— >  X
         rect:     RECT of floats                   bottom | ┌─┐
        this value shows edges of area                 top | └─┘RECT
                                                        Y  V             // carthesian system is computer like (OY goes down)  
         conds:    RECT of lambdas                       
        conditions on the rect edges

         x_step:   float
        size of grid on Ox

         y_step:   float
        size of grid on Oy
        if was not set then = x_step

         accuracy: float
        accuracy of lin alg calculations

    Return:

         array
        the array with counted values on points of our grid inside the rect
        they are numbered from 0 to i_max and does upside down according to the picture, from left to right

    """         


    # check edges and steps    TODO: mb delete
    assert int((rect.right - rect.left) / x_step) == ((rect.right - rect.left) / x_step),  "incorrect x_step"
    assert int((rect.top - rect.bottom) / y_step) == ((rect.top - rect.bottom) / y_step),  "incorrect y_step"

    # calculate max values for j, k and i
    n = int( (rect.right - rect.left) / x_step ) + 1
    m = int( (rect.top - rect.bottom) / y_step ) + 1
    i_max = (n - 2) * (m - 2)
    if (DEBUG):
        print("\n\nBoundValPuassonRectEq_DEBUG\nn =", n, " m =", m, " i_max =", i_max, "\nx_step = ", x_step, "\ny_step = ", y_step, end="\n")

    # macros to calculate one coords from another
    j2x  =  lambda j     :  j * x_step + rect.left
    k2y  =  lambda k     :  k * y_step + rect.bottom
    kj2i =  lambda k, j  :  (k - 1) * (n - 2) + j - 1
    i2jk =  lambda i     :  ( i % (n - 2) + 1,  int(i / (n - 2)) + 1 )

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
            rvec[i] -= conds.right([j2x(j + 1), k2y(k)]) / x_step**2
        else:
            matr[i][kj2i(k, j+1)] = 1 / x_step**2

        # u_j-1_k
        if (j - 1 == 0):
            rvec[i] -= conds.left([j2x(j - 1), k2y(k)]) / x_step**2
        else:
            matr[i][kj2i(k, j-1)] = 1 / x_step**2

        # u_j_k+1
        if (k + 1 == m - 1):
            rvec[i] -= conds.top([j2x(j), k2y(k + 1)]) / y_step**2
        else:
            matr[i][kj2i(k+1, j)] = 1 / y_step**2

        # u_j_k-1
        if (k - 1 == 0):
            rvec[i] -= conds.bottom([j2x(j), k2y(k - 1)]) / y_step**2
        else:
            matr[i][kj2i(k-1, j)] = 1 / y_step**2


    # print the matrix
    if (DEBUG):
        _debug_str = ""
        _debug_str += "\nThe matrix of linear equation system:\n\n"
        table_col_size = GetIntDigitsCount(matr[0][0]) + 3
        _debug_str += " "*table_col_size
        for j in range(i_max):
            _debug_str += str(j + 1).center(table_col_size)
        _debug_str += "\n"
        for k in range(i_max):
            _debug_str += str(k + 1).center(table_col_size)
            for j in range(i_max):
                if(matr[k][j]):
                    _debug_str += str(round(matr[k][j], 2)).rjust(table_col_size)
                else:
                    _debug_str += " "*table_col_size
            _debug_str += " |" + str(round(rvec[k], 2)).rjust(table_col_size) + "\n"
        if (i_max < 20):
            print(_debug_str)
        else:
            print("Lin System matrix u can see in the file \"eq sys matr.txt\"\n")
        f = open("eq sys matr.txt", "w")
        try:
            f.write(_debug_str)
        except BaseException:
            TextException()
            print("FILE DEBUG ERROR")
        finally:
            f.close()
        print()


    """
    Solve the system
    """
    x = np.array([1 for i in range(i_max)])
    x = LinAlgRelax(x, matr, rvec, accuracy, 0.9)

    return x
