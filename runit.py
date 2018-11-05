from random import randint
import gauss
import algoritmos as al
import numpy as np

# T = np.arange(0, 30)
T = np.arange(0, 10)
h = al.h(T)
mi = al.mi_array(h)
lambda_array = al.lambda_array(h)
matrix = al.make_matrix(len(h), mi, lambda_array)


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
