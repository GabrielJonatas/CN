import calc
import solver
import numpy as np

T = calc.T

x00 = calc.x00()
xnn = calc.xnn(28)

h = calc.h(T)
n = len(h)

# Vetores mi e Lambda & matriz do sistema (2)
mi = calc.mi_array(h)
lambda_array = calc.lambda_array(h)
matrix = calc.make_matrix(n, mi, lambda_array)

# Resolve por Gauss
M_values = calc.M_gauss(n, x00, xnn, mi, lambda_array, h)

M = np.zeros(len(M_values))
i = 0
for k in M_values:
    M[i] = float(M_values[k])
    i += 1


B = calc.B(M, h, n)
A = calc.A(M, h, n)

d = calc.make_d(h, x00, xnn)


print('\n========================= TAREFA 1 ATIVIDADE 1 =========================')
estritamente_dominante = calc.is_strictly_diagonal_dominant(matrix)
print(f'É estritamente dominante? --> {estritamente_dominante}')
sigma = calc.make_sigma(matrix)
print(f'Valor de Sigma: {sigma}')

print('\n========================= TAREFA 1 | ATIVIDADES 2 & 4 =========================')
print('OBS: Calculamos a estimativa de número mínimo de iterações para a tolerância desejada DENTRO do método de Jacobi.')
# Resolve por Jacobi
jacobi = solver.jacobi(matrix, d, mi, lambda_array)
print(f'\nSOLUCAO OBTIDA POR JACOBI:\n{jacobi}')
print(f'\nSOLUCAO OBTIDA POR GAUSS:\n{M}')

print('\n========================= TAREFA 1 ATIVIDADE 3 =========================')
norm = solver.infinity_norm(jacobi, M)
print(f'\nNORMA INFINITA --> || JACOBI - GAUSS ||∞ = {norm}')

print('\n========================= TAREFA 2 ATIVIDADE 1 =========================')
f205 = calc.f(20.5, h, M, A, B, 13)
f21 = calc.f(21, h, M, A, B, 13)
print(f'Intervalo [20.5, 21]')
print(f'f(20.5)*f(21) = {f205 * f21} > 0')
calc.newton(21, 10**-8, A, B, M, 13, h)
