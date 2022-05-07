"""
The sample solution of Boundary value problem with Puasson equation
"""

from params import *
from base_objs import *
from solution_methods import *

import numpy as np
import scipy as sp
import scipy.linalg as lalg

import matplotlib.pyplot as plot
import matplotlib.patches as patches

import time



def main():

    """
Solve the problem using our method
    """

    start_time = time.time()
    x = BoundValPuassonRectEq(f, rect, conds, h, l)
    end_time   = time.time()

    print("Calculation was done in", round(end_time - start_time, 3), "seconds\n")

    """
Exact solution
    """
    n = int( (rect.right - rect.left) / h ) + 1
    i_max = (n - 2) ** 2
    ex_x = np.zeros(i_max)
    for i in range(i_max):
        j, k = ( i % (n - 2) + 1,  int(i / (n - 2)) + 1 )
        ex_x[i] = ex_sol([j * h, k * l])

    print("\n\n\
Lets look at the solutions\
    ")
    table_col_width = 11
    print("┌────┬─", end="")
    print("─"*table_col_width, "─"*table_col_width, sep="─┬─", end="─┐\n")
    print(
        "│" + "i".rjust(3),
        "Low Relax".center(table_col_width),
        "Exact".center(table_col_width),
        sep=" │ ", end=" │\n"
    )
    for i in range(i_max):
        print("├────┼─", end="")
        print("─"*table_col_width, "─"*table_col_width, sep="─┼─", end="─┤\n")
        print(
            "│" + str(i).rjust(3),
            str(round(x[i],       table_col_width - 3)).center(table_col_width),
            str(round(ex_x[i],    table_col_width - 3)).center(table_col_width),
            sep=" │ ", end=" │\n"
        )
    print("└────┴─", end="")
    print("─"*table_col_width, "─"*table_col_width, sep="─┴─", end="─┘\n")



    """
And now we will get matrix of so m matrix which hmm we need the matrix ...    ok?
    """

    final = np.zeros((n, n))

    # fill edge conditions

    for j in range(n):
        final[0][j]     = conds.left([j * h, rect.bottom])
    for k in range(n):
        final[k][0]     = conds.bottom([rect.left, k * l])
    for j in range(n):
        final[n - 1][j] = conds.top([j * h, rect.top])
    for k in range(n):
        final[k][n - 1] = conds.right([rect.right, k * l])

    # fill our solution

    for i in range(i_max):
        j, k = ( i % (n - 2) + 1,  int(i / (n - 2)) + 1 )
        final[k][j] = x[i]

    print("Solved grid of temperatures")
    print(final)


    """
    Lets see in graph
    """

    try:

        # for colours
        temp_max   = final.max()
        temp_min   = final.min()
        temp_delta = temp_max - temp_min

        fig, ax = plot.subplots()

        for k in range(n):
            for j in range(n):

                theColor = (0,0,0, (temp_max - final[k][j]) / temp_delta * 0.9 + 0.1)

                ax.add_patch(
                    patches.Rectangle(
                        (j * h, k * l),
                        h,
                        l,
                        facecolor = theColor,
                        fill=True
                    ) )

        plot.show()

    except BaseException:
        print("\n!!! PLOT ERROR\n")
        print(TextException())




if __name__ == "__main__":
    main()
    print("-"*200, "мучас грасиас офишион ешто паравасотрощ шууууууууууу".upper(), "="*200, sep="\n", end="\n"*13)
