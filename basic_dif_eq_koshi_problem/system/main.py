# solution run

# import params
from params import *
# import process func
from solution_methods import *

import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt


DISABLE_PLOTS            = FALSE
TURN_ON_THE_COLORS       = TRUE
PLOT_DEBUG_POINTS        = FALSE
PLOT_GRID_XLINESSS_COUNT = 20
PLOT_MARKER_SIZE         = 6

def ShowSolutionsPlots(x: list, title: str, labels: list = False):

    """
    Show plot of your solutions

    x is list of many solutions where each one is list of dots (x_2,y,z) (list of list of list of float)
    """

    if (DISABLE_PLOTS):
        return


    # colors for different plots (if TURN_ON_THE_COLORS is TRUE)
    colors_for_plot = [
        "red",
        "blue",
        "green",
        "purple"
    ]



    # restructure vectors of dots from [ [x1, y1], ... ] to [x1, ...], [y1, ...]

    x_new = []

    for i in range(len(x)):

        x_step = [[], []]

        for k in range(len(x[i])):
            x_step[0].append(x[i][k][0])
            x_step[1].append(x[i][k][1])

        x_new.append(x_step)



    try:

        fig, ax = plt.subplots()    

        fig.set_size_inches(13, 9) 

        range_start = x[0][0][0]
        range_end   = x[0][len(x[0])-1][0]
        grid_sizze  = (range_end - range_start) / PLOT_GRID_XLINESSS_COUNT
        ax.set_xticks(np.arange( range_start,  range_end + grid_sizze, grid_sizze))
        ax.grid(color='black', linewidth=0.5, linestyle='--')


        for i in range(len(x)):

            if (TURN_ON_THE_COLORS):
                try:
                    theColor = colors_for_plot[i]
                except BaseException:
                    theColor = (0,0,0, 1 / len(x) * (i+1) )
            else:
                theColor = (0,0,0, 1 / len(x) * (i+1) )

            line, = ax.plot(x_new[i][0], x_new[i][1], color=theColor, marker='o', markersize=PLOT_MARKER_SIZE)
            if (labels):
                line.set_label(labels[i])


        plt.title(title)
        fig.canvas.manager.set_window_title("Solution of differetintial equations' system")
        ax.legend(labels)

        plt.show()

    except BaseException:
        print("\n!!! PLOT ERROR\n")
        print(TextException())



# here the calculations start
def main():

    # check if data is correct
    assert RANGE.end > RANGE.start,   "Incorrect RANGE (end < start)"
    assert len(FUNC) == len(X0) - 1,  "Incorrect problem vectors (count of function != size of Koshi task)"


    step = (RANGE.end - RANGE.start) / SEGMENTS_COUNT



    """
    Looking for solution using 2-order Runge-Kutta method
    """

    try:

        x_2  = RungeKutta2(FUNC, X0, step,   SEGMENTS_COUNT)
        x2_2 = RungeKutta2(FUNC, X0, step/2, SEGMENTS_COUNT*2)
        x4_2 = RungeKutta2(FUNC, X0, step/4, SEGMENTS_COUNT*4)
        x8_2 = RungeKutta2(FUNC, X0, step/8, SEGMENTS_COUNT*8)
        ShowSolutionsPlots([x_2, x2_2, x4_2, x8_2], "Runge-Kutta 2-order method", ["full step", "half step", "quarter step", "1/8 step"])

    except BaseException:
        print(TextException())



    """
    Now looking for solution using 3-order Runge-Kutta method
    """

    try:

        x_3  = RungeKutta3(FUNC, X0, step,   SEGMENTS_COUNT)
        x2_3 = RungeKutta3(FUNC, X0, step/2, SEGMENTS_COUNT*2)
        x4_3 = RungeKutta3(FUNC, X0, step/4, SEGMENTS_COUNT*4)
        x8_3 = RungeKutta3(FUNC, X0, step/8, SEGMENTS_COUNT*8)
        ShowSolutionsPlots([x_3, x2_3, x4_3, x8_3], "Runge-Kutta 3-order method", ["full step", "half step", "quarter step", "1/8 step"])

    except BaseException:
        print(TextException())



    """
    Now looking for solution using 4-order Runge-Kutta method
    """

    try:

        x_4  = RungeKutta4(FUNC, X0, step,   SEGMENTS_COUNT)
        x2_4 = RungeKutta4(FUNC, X0, step/2, SEGMENTS_COUNT*2)
        x4_4 = RungeKutta4(FUNC, X0, step/4, SEGMENTS_COUNT*4)
        x8_4 = RungeKutta4(FUNC, X0, step/8, SEGMENTS_COUNT*8)
        ShowSolutionsPlots([x_4, x2_4, x4_4, x8_4], "Runge-Kutta 4-order method", ["full step", "half step", "quarter step", "1/8 step"])

    except BaseException:
        print(TextException())



    # SciPy solution
    try:
        
        print("\n\nSolving using SciPy")
        def _f4sp(y: list, x: float) -> list:
            return [
                FUNC[0]([ x, y[0], y[1] ]),
                FUNC[1]([ x, y[0], y[1] ]),
            ]
        _x4sp = np.linspace(RANGE.start, RANGE.end, SEGMENTS_COUNT+1)
        _x_sp = odeint(_f4sp, X0[1:], _x4sp)

        x_sp = []
        for i in range(len(_x4sp)):
            x_sp.append([ _x4sp[i], _x_sp[i][0], _x_sp[i][1] ])

    except BaseException:
        print(TextException())



    print("\n\n\
    And here we can see da difference beetwin the 2,3 and 4 -order meth(ods)                    \n\
    with only 11 steps solution so we can see it                                                \n\
    and also the scipy solution                                                                 \n\
    ")

    ShowSolutionsPlots([x_2, x_3, x_4, x_sp], "Runge-Kutta 2+3+4-order methods (11 steps) and SciPy", ["2-order method", "3-order method", "4-order method", "math library solution"])




if __name__ == "__main__":
    main()
    print("-"*200, "мучас грасиас офишион ешто паравасотрощ шууууууууууу".upper(), "="*200, sep="\n", end="\n"*13)
