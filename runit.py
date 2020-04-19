import algoritmos
import solver
import numpy as np

T = algoritmos.T

x00 = algoritmos.x00()
xnn = algoritmos.xnn(28)

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
print(f'É estritamente dominante? --> {estritamente_dominante}')
sigma = algoritmos.make_sigma(matrix)
print(f'Valor de Sigma: {sigma}')

print('\n========================= TAREFA 1 | ATIVIDADES 2 & 4 =========================')
print('OBS: Calculamos a estimativa de número mínimo de iterações para a tolerância desejada DENTRO do método de Jacobi.')
# Resolve por Jacobi
jacobi = solver.jacobi(matrix, y, mi, lambda_array)
print(f'\nSOLUCAO OBTIDA POR JACOBI:\n{jacobi}')
print(f'\nSOLUCAO OBTIDA POR GAUSS:\n{M}')

print('\n========================= TAREFA 1 ATIVIDADE 3 =========================')
norm = solver.infinity_norm(jacobi, M)
print(f'\nNORMA INFINITA --> || JACOBI - GAUSS ||∞ = {norm}')

print('\n========================= TAREFA 2 ATIVIDADE 1 =========================')
f205 = algoritmos.f(20.5, h, M, A, B, 13)
f21 = algoritmos.f(21, h, M, A, B, 13)
print(f'Intervalo [20.5, 21]')
print(f'f(20.5)*f(21) = {f205 * f21} > 0')
algoritmos.newton(21, 10**-8, A, B, M, 13, h)
