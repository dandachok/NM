import numpy as n
from math import *
s = 3
eps = 0.3
A = n.array([[9,2,-7],
[2,-4,-1],
[-7,-1,1]], float)
R = n.eye(3)
print("Cобственные значения python", n.linalg.eig(A)[0])
print("Cобственные векторы python\n", n.linalg.eig(A)[1])
while 1:
    max_dig = abs(A[0][1])
    max_i = 0
    max_j = 1
    for i in range(s):
        for j in range(s):
            if max_dig < abs(A[i][j]) and i != j:
                max_dig = abs(A[i][j])
                max_i = i
                max_j = j
    U = n.eye(s)
    fi = 0.5*atan((2*max_dig)/(A[max_i][max_i]-A[max_j][max_j]))
    U[max_i][max_j] = -sin(fi)
    U[max_i][max_i] = U[max_j][max_j] = cos(fi)
    U[max_j][max_i] = sin(fi)
    R = n.dot(R,U)
    print(U)
    A = n.dot(n.dot(n.transpose(U), A),U)
    t = 0
    for l in range(s):
        for m in range(l):
            t = t + pow(A[l][m],2)
    t = sqrt(t)
    if t < eps:
        break