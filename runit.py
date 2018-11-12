import algoritmos
import solver
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

y = algoritmos.make_d(h, x00, xnn)


print('\n========================= TAREFA 1 ATIVIDADE 1 =========================')
estritamente_dominante = algoritmos.is_strictly_diagonal_dominant(matrix)
print('É estritamente dominante? --> {}'.format(estritamente_dominante))
sigma = algoritmos.make_sigma(matrix)
print('Valor de Sigma: {}'.format(sigma))

print('\n========================= TAREFA 1 | ATIVIDADES 2 & 4 =========================')
print('OBS: Calculamos a estimativa de número mínimo de iterações para a tolerância desejada DENTRO do método de Jacobi.')
# Resolve por Jacobi
jacobi = solver.jacobi(matrix, y, mi, lambda_array)
print('\nSOLUCAO OBTIDA POR JACOBI:\n{}'.format(jacobi))
print('\nSOLUCAO OBTIDA POR GAUSS:\n{}'.format(M))

print('\n========================= TAREFA 1 ATIVIDADE 3 =========================')
norm = solver.infinity_norm(jacobi, M)
print('\nNORMA INFINITA --> || JACOBI - GAUSS ||∞ = {}'.format(norm))

print('\n========================= TAREFA 2 ATIVIDADE 1 =========================')
result = newton(u_0, epsilon, A, B, M, 13, h)