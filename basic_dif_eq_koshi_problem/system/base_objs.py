# basic objects for calculations

from math import fabs
from xml.etree.ElementTree import tostring
import numpy as np
import matplotlib.pyplot as plot

import sys
import linecache


TRUE  = 1
FALSE = 0


class Range:

    def __init__(self):
        self.start  = 0
        self.end    = 0

    def __init__(self, start: float, end: float):
        self.start  = start
        self.end    = end

    def to_str(self) -> str:
        return "(" + str(self.start) + ", " + str(self.end) + ")"

    def contains(self, x: float) -> bool:
        return (x >= self.start) and (x <= self.end)



def TextException():

    """
    exception text compiler
    """
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    return 'EXCEPTION IN ({}, LINE {} "{}"):\n{}'.format(filename, lineno, line.strip(), exc_obj)



# Accuracy offset
ACCURACY_ROUND_OFFSET = 5
# u can increase count of symbols after point using this constant

def GetIntDigitsCount(x: float) -> int:
    """
    Returns count of digits before point
    """
    return len(str(x).split('.')[0]) + 1#for sign

def GetFloatDigitsCount(x: float) -> int:
    """
    Returns count of digits after point
    """
    try:
        return len(str(x).split('.')[1])
    except BaseException:
        try:
            return int(str(x).split('-')[1])
        except BaseException:
            return 0


def RoundBySymbCount(x: float, symb_count: float) -> float:
    """
    Rounds the float X according to given symbols after point count
    """
    return float(f"{x:.{symb_count}f}")    

def RoundByAccur(x: float, accuracy: float) -> float:
    """
    Rounds the float X according to given accuracy
    """
    return RoundBySymbCount(x, GetFloatDigitsCount(accuracy) + ACCURACY_ROUND_OFFSET)

def RoundByAccurVec(x: list, accuracy: float) -> list:
    """
    Rounds the float X according to given accuracy
    """
    t = []
    for i in x:
        t.append(float(f"{i:.{GetFloatDigitsCount(accuracy)}f}"))
    return t



def CalculateVectorFunction(f: list, x: list) -> list:
    """
    Calculates value of vector function in some vector point
    """
    f_value = []
    for i in f:
        f_value.append(i(x))
    return f_value

def VectorMax(v: list):
    max = v[0]
    for i in v:
        if i > max:
            max = i
    return max

def VectorsMaxDelta(v1: list, v2: list):
    """
    Calculates deltas by every axis and returns the max of them
    """
    if len(v1) != len(v2):
        raise "Incorrect arguements"

    deltas = []
    for i in range(len(v2)):
        deltas.append(fabs(v2[i] - v1[i]))
    
    return VectorMax(deltas)


def ShowRootPlot(f, root: float, range: Range, accurasy: float):

    """
    Show plot of your func on the arry with root point

    Use only numpy math functions
    """

    try:
        t = np.arange(range.start, range.end, accurasy)
        fig, ax = plot.subplots()
        ax.set_title("Root: " + str(root) + " +- " + str(accurasy))
        plot.grid(which='minor', color = 'k', linestyle = ':')
        plot.grid(which='major', color = 'k', linewidth = 1,linestyle = ':')
        plot.plot(t, f(t), "r--")
        plot.plot([root], [f(root)], "bo")
        plot.show()
    except BaseException:
        print("\n!!! PLOT ERROR\n")
        print(TextException())
