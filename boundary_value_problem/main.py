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
import matplotlib.ticker as ticker


import ctypes
lib = ctypes.windll.kernel32



def main():

    """
Solve the problem using our method
    """
    
    print("Calculations started...")

    start_time = lib.GetTickCount64()
    x = BoundValPuassonRectEq(f, rect, conds, h, l, accuracy = accuracy)
    end_time   = lib.GetTickCount64()

    print("\n\nCalculations were done in", (end_time - start_time) / 1000, "seconds\n")

    """
Exact solution
    """
    n = int( (rect.right - rect.left) / h ) + 1
    m = int( (rect.top - rect.bottom) / l ) + 1
    i_max = (n - 2) * (m - 2)
    ex_x = np.zeros(i_max)
    for i in range(i_max):
        j, k = ( i % (n - 2) + 1,  int(i / (n - 2)) + 1 )
        ex_x[i] = ex_sol([j * h + rect.left, k * l + rect.bottom])

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
    index = 0
    while (index < i_max):
        print("├────┼─", end="")
        print("─"*table_col_width, "─"*table_col_width, sep="─┼─", end="─┤\n")
        print(
            "│" + str(index + 1).rjust(3),
            str(round(x[index],       table_col_width - 3)).center(table_col_width),
            str(round(ex_x[index],    table_col_width - 3)).center(table_col_width),
            sep=" │ ", end=" │\n"
        )
        index += 1
        if (i_max > 50 and index == 26):
            print("\n    ...\n")
            index = i_max - 25
    print("└────┴─", end="")
    print("─"*table_col_width, "─"*table_col_width, sep="─┴─", end="─┘\n")

    print("Final accuracy:", VectorsMaxDelta(x, ex_x), "\n")



    """
And now we will get matrix of so m matrix which hmm we need the matrix ...    ok?
    """

    final = np.zeros((m, n))

    # fill edge conditions

    for j in range(n):
        final[0][j]     = conds.bottom([j * h + rect.left, rect.bottom])
    for k in range(m):
        final[k][0]     = conds.left([rect.left, k * l + rect.bottom])
    for j in range(n):
        final[m - 1][j] = conds.top([j * h + rect.left, rect.top])
    for k in range(m):
        final[k][n - 1] = conds.right([rect.right, k * l + rect.bottom])

    # fill our solution

    for i in range(i_max):
        j, k = ( i % (n - 2) + 1,  int(i / (n - 2)) + 1 )
        final[k][j] = x[i]

    print("\nSolved grid of temperatures")

    table_col_width = 6
    table_col_count = n
    print("┌─", ("─"*table_col_width + "─┬─")*(table_col_count-1), "─"*table_col_width + "─┐", sep="")
    for k in range(m):
        if (k > 0):
            print("├─", ("─"*table_col_width + "─┼─")*(table_col_count-1), "─"*table_col_width + "─┤", sep="")
        print("│", end="")
        for el in final[k]:
            print(str(round(el, table_col_width - 2)).center(table_col_width + 2), "│", sep="", end = "")
        print()
    print("└─", ("─"*table_col_width + "─┴─")*(table_col_count-1), "─"*table_col_width + "─┘", sep="")




    """
    Lets see in graph
    """

    try:

        # for colours
        temp_max   = final.max()
        temp_min   = final.min()
        temp_delta = temp_max - temp_min

        fig, ax = plot.subplots()
        ax.set_xlim([rect.left, rect.left + n*h])
        ax.set_ylim([rect.bottom, rect.bottom + m*l])

        for k in range(m):
            for j in range(n):

                theColor = (0,0,0, (temp_max - final[k][j]) / temp_delta * 0.9 + 0.1)

                ax.add_patch(
                    patches.Rectangle(
                        (j * h + rect.left, k * l + rect.bottom),
                        h, l,
                        facecolor = theColor,
                        fill=True
                    ))

        plot.show()

    except BaseException:
        print("\n!!! PLOT ERROR\n")
        print(TextException())




if __name__ == "__main__":
    main()
    print("-"*200, "мучас грасиас офишион ешто паравасотрощ шууууууууууу".upper(), "="*200, sep="\n", end="\n"*13)
