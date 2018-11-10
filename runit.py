import algoritmos
import numpy as np

T = algoritmos.T

x00 = algoritmos.x00(T)
xnn = algoritmos.xnn(T, 28)

h = algoritmos.h(T)
n = len(h)

# Vetores mi e Lambda & matriz do sistema (2)
mi = algoritmos.mi_array(h)
lambda_array = algoritmos.lambda_array(h)
matrix = algoritmos.make_matrix(n, mi, lambda_array)

# Resolve por Gauss
M_values = algoritmos.M_gauss(n, x00, xnn, mi, lambda_array, h)

M = np.zeros(len(M_values))
i = 0
for k in M_values:
    M[i] = float(M_values[k])
    i += 1

B = algoritmos.B(M, h, n)
A = algoritmos.A(M, h, n)

