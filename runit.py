import algoritmos as al
import numpy as np

x00 = al.x00
xnn = al.xnn(28)

h = al.h(al.T)
n = len(h)

mi = al.mi_array(h)
lambda_array = al.lambda_array(h)
matrix = al.make_matrix(n, mi, lambda_array)

M_values = al.M_gauss(n, x00, xnn, mi, lambda_array, h)

M = np.zeros(len(M_values))
i = 0
for k in M_values:
    M[i] = float(M_values[k])
    i += 1

B = al.B(M, h, n)
A = al.A(M, h, n)

