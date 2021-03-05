import numpy as np
import math as m
eps = 10**(-5)
def f_f(x1, x2):
    return (x1**2 + 16)*x2 - 64
    #return x1 - m.cos(x2) - 1
    #return 0.1* x1**2 + x1 + 0.2* x2**2 - 0.3

def s_f(x1, x2):
    return (x1 - 2)**2 + (x2 - 2)**2 - 16
    #return x2 - m.log(x1 + 1) - 3
    #return 0.2* x1**2 + x2 - 0.1*x1*x2 - 0.7

def f_f_f_x1(x1,x2):
    return 2*x1*x2
    #return 1
#return 0.2* x1 + 1

def f_f_f_x2(x1,x2):
    return x1**2 + 16
    #return m.sin(x2)
#return 0.4* x2

def s_f_f_x1(x1,x2):
    return 2*x1 - 4
    #return -(1/ (x1 + 1))
#return 0.4*x1 - 0.1*x2

def s_f_f_x2(x1,x2):
    return 2*x2 - 4
    #return 1
#return 1 - 0.1*x1

def J(x1,x2):
    return np.linalg.det([[f_f_f_x1(x1,x2), f_f_f_x2(x1,x2)],[s_f_f_x1(x1,x2), s_f_f_x2(x1,x2)]] )

def A1(x1,x2):
    return np.linalg.det([[f_f(x1,x2), f_f_f_x2(x1,x2)],[s_f(x1,x2), s_f_f_x2(x1,x2)]])

def A2(x1,x2):
    return np.linalg.det([[f_f_f_x1(x1,x2), f_f(x1,x2)],[s_f_f_x1(x1,x2), s_f(x1,x2)]])

def norm(x1,x2, x1_prev, x2_prev):
    return abs(abs(x1) - abs(x1_prev)) + abs(abs(x2) - abs(x2_prev))


def Newton():
    x1_prev = 5
    x2_prev = 1
    x1 = 5.5
    x2 = 1.1
    while norm(x1, x2, x1_prev, x2_prev) > eps:
        x1_prev = x1
        x2_prev = x2
        detA1 = A1(x1,x2).copy()
        detA2 = A2(x1,x2).copy()
        detJ = J(x1,x2).copy()
        x1 = x1 - detA1/detJ
        x2 = x2 - detA1/detJ
    print("Ответ:", x1, x2)