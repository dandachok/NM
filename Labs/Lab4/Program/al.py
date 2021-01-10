import numpy as np
import math as m
import matplotlib.pyplot as plt
left = 0.
right = 1.
h = 0.1
yo = 1.
zo = 1.
flag = 0
z_for_Adams = []
def f_real(x, y, z):
    return (1 + x) / m.exp(x**2)
def f(x, y, z):
    return z
def g(x, y, z):
    return -(4*x**2 + 2)*y - 4*x*z
def count_RR(yh, y2h):
    return yh+(yh-y2h)/(2**4 - 1)
def RR(xi, y1, y2, h):
#у у2 шаг меньше, а xi это массив значений для шага h/2, то есть это массив x для y1.
    ans = []
    i = 0
    while i < len(xi):
        ans.append(count_RR(y1[i],y2[i//2]))
        i += 2
    return ans

def count_k(h, x, y, z):
    return h*f(x, y, z)

def count_l(h, x, y, z):
    return h*g(x,y,z)

def grow_y(K):
    return (1/6)*(K[-1] + K[0] + 2*sum(K[1:(len(K)-1)]))
def grow_z(L):
    return (1/6)*(L[-1] + L[0] + 2*sum(L[1:(len(L)-1)]))

def RK(h):
    xi = np.arange(left, right + h, h)
    yi = []
    yi.append(yo)
    zi = []
    zi.append(zo)
    yi_real = []
    for i in range(len(xi)):
        K = []
        L = []
        for j in range(4):
            if j == 0:
                K.append(count_k(h, xi[i], yi[i], zi[i]))
                L.append(count_l(h, xi[i], yi[i], zi[i]))
            elif j == 3:
                K.append(count_k(h, xi[i] + h, yi[i] + K[j-1], zi[i] + L[j-1]))
                L.append(count_l(h, xi[i] + h, yi[i] + K[j-1], zi[i] + L[j-1]))
            else:
                K.append(count_k(h, xi[i] + 0.5*h, yi[i] + 0.5*K[j-1], zi[i] + 0.5*L[j-1]))
                L.append(count_l(h, xi[i] + 0.5*h, yi[i] + 0.5*K[j-1], zi[i] + 0.5*L[j-1]))
        print("K: ", K)
        print("y: ", yi)
        dyi = grow_y(K)
        dzi = grow_z(L)
        yi_next = yi[i] + dyi
        zi_next = zi[i] + dzi
        yi.append(yi_next)
        zi.append(zi_next)
        yi_real.append(f_real(xi[i], yi[i], zi[i]))

    yi = yi[:len(xi)]
    return yi, zi