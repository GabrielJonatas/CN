from random import randint
import gauss
import algoritmos as al
import numpy as np

x00 = al.x00
xnn = al.xnn(28)

h = al.h(al.T)
n = len(h)

mi = al.mi_array(h)
lambda_array = al.lambda_array(h)
matrix = al.make_matrix(n, mi, lambda_array)
d = al.d(h, x00)

M = al.M_gauss(n, x00, xnn, mi, lambda_array)

B = al.B(M, h, n)
A = al.A(M, h, n)


shape = (3, 3)
A = np.zeros(shape)

for i in range(shape[0]):
    for j in range(shape[1]):
        A[i][j] = randint(1, 9)

b = np.zeros(shape[0])
for i in range(shape[0]):
    b[i] = randint(1, 9)

result = gauss.for_matrix(A)
print(result)

print('----------------------')
b = list(b)
A = list(A)
gauss.linearsolver(A, b)
