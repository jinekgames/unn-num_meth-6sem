"""
Entry params for sample solution
"""


from base_objs import *
import numpy as np


# rect edges            a b    
_A = 0      #         +——————— >  X
_B = 1      #       c | ┌─┐
_C = 0      #       d | └─┘
_D = 1      #      Y  V        // carthesian system is computer like (OY goes down)

rect = RECT_DOUBLE(_A, _D, _B, _C)


# x and y step
h = 1/4
l = 1/4


# edge conditions
x_y0 =   lambda x :   1                  #   u = u(x, c)            
x1_y =   lambda x :   np.exp( -x[1] )    #   u = u(b, y)         
x_y1 =   lambda x :   np.exp( -x[0] )    #   u = u(x, d)     
x0_y =   lambda x :   1                  #   u = u(a, y)

conds = RECT_FUNC(x0_y, x_y1, x1_y, x_y0)


# function from equations (right part of the Puasson equation)
f    =   lambda x :   np.exp( -x[0]*x[1] ) * (x[0]**2 + x[1]**2)


# exact solution (the u = u(x, y) function we want to get)
def ex_sol(x: list) -> list:
    return np.exp(-x[0]*x[1])
