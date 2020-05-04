import numpy as np
# from sympy import Matrix, linsolve, solve_linear_system, symbols

# Dados da queda livre
tao = 23  # Numero do grupo, usado na busca binaria
distance = 100 - 2.3 * tao
T = np.arange(0, 30)  # Vetor de tempo em segundos
X = [
    0, 1, 2.4, 4.1, 6, 8.2, 10.6, 13.4, 16.4, 19.7,
    23.3, 27, 31.2, 35.5, 40.1, 45, 50.2, 55.6, 61.3, 67.3,
    73.6, 80.1, 86.9, 94, 101.3, 109, 116.9, 125, 133.4, 142.1
]  # Distância percorrida em queda livre ao longo do tempo

epsilon = 10**-3


def x00():
    return ((X[2] - X[1]) / (T[2] - T[1]) - (X[1] - X[0]) / (T[1] - T[0])) / (T[2] - T[0])


def xnn():
    n = len(T) - 1
    return ((X[n] - X[n - 1]) / (T[n] - T[n - 1]) - (X[n - 1] - X[n - 2]) / (T[n - 1] - T[n - 2])) / (T[n] - T[n - 2])


# Tarefa 1, Exercicio 1
# Determinação de dominancia diagonal estrita de uma matriz
def is_strictly_diagonal_dominant(A):
    for i in range(len(A)):
        line_sum = 0

        for j in range(len(A[0])):
            if j != i:
                line_sum += A[i][j]

        if A[i][i] < line_sum:
            print(f'Linha {i} NÃO é diagonal dominante.\t\tA[{i}][{i}] = {A[i][i]}\tSoma da linha = {line_sum}')
            return False

    # Passou por todas linhas e nenhuma soma deu maior que a diagonal
    return True


# Tarefa 1, Exercicio 1
# Calculo de sigma
def make_sigma(A):
    sigma = np.zeros(len(A))

    for i in range(len(A)):
        sigma[i] = (1 / abs(A[i][i])) * (sum(A[i]) - A[i][i])

    return max(sigma)


def h(T):
    n = len(T)
    h = np.zeros(n - 1)

    # Valores de h são as diferenças entre valores de t
    for i in range(n - 1):
        h[i] = T[i + 1] - T[i]
    return h


# VALORES DE µ (mi)
def mi_array(h):
    n = len(h) - 1
    mi = np.zeros(n)

    mi[0] = 0

    for i in range(1, n):
        mi[i - 1] = h[i] / (h[i] + h[i + 1])

    return mi


def lambda_array(h):
    n = len(h) - 1
    lambda_array = np.zeros(n)

    lambda_array[0] = 0

    for i in range(1, n):
        lambda_array[i] = h[i + 1] / (h[i] + h[i + 1])

    return lambda_array


def make_matrix(n, mi_array, lambda_array):
    # Instancia matrix nxn
    A = np.zeros((n, n))

    for i in range(n):
        A[i][i] = 2
        if i != (n - 1):
            A[i][i + 1] = lambda_array[i]
            A[i + 1][i] = mi_array[i]

    return A


def make_d(h, x00, xnn):
    n = len(h)
    d = np.zeros(n)

    d[0] = 2 * x00
    d[n - 1] = 2 * xnn

    for i in range(1, n - 1):
        d[i] = (6 / (h[i] + h[i + 1])) * ((X[i + 1] - X[i] / h[i + 1]) - (X[i] - X[i - 1]) / h[i])

    return d


def B(M, h, n):
    B = np.zeros(n - 1)

    for i in range(n - 1):
        B[i] = X[i] - (M[i] / 6) * h[i + 1]**2

    return B


def A(M, h, n):
    A = np.zeros(n - 1)

    for i in range(n - 1):
        A[i] = (X[i + 1] - X[i]) / h[i + 1] - ((M[i + 1] - M[i]) / 6) * h[i + 1]

    return A


# VETOR DE VALORES M
def M_gauss(n, x00, xnn, mi, lambda_array, h):
    A = make_matrix(n, mi, lambda_array)
    d = make_d(h, x00, xnn)
    return gauss(A, d)


# IMPLEMENTA MÉTODO DE GAUSS
def gauss(A, y, m=0):
    if m == 0:
        m = len(A)
    for j in range(1, m):  # para j de 1 até m-1
        find_pivot(A, y, j, m)
    # print("\n")
    # print(A)
    for i in range(j+1, m+1):  # para i de 1 até m
        mi = - (element(A, i, j) / element(A, j, j))
        for k in range(j, m+1): # para k de j até m
            set_element(A, i, k, element(A, i, k) + mi * element(A, j, k))
        y[i-1] = y[i-1] + mi * y[j-1]
        # print("\n")
        # print(A)
    x = np.zeros(len(y))  # Cria array e preenche com zeros
    i = m
    while (i >= 1):
        x[i-1] = y[i-1]
        for j in range(i+1, m+1):  # para j de i+1 até m
            x[i-1] = x[i-1] - element(A, i, j) * x[j-1]
        x[i-1] = x[i-1] / element(A, i, i)
        i -= 1
    # print(x)
    return x

def find_pivot(A, y, j, m):
    if element(A, j, j) == 0:
        matrix_singular = True
        for k in range(j+1, m+1): # para k de j+1 até m
            if element(A, k, j) != 0:
                matrix_singular = False
                swap_lines(A, j, k)
                swap_lines(y, j, k)
                break
        if matrix_singular:
            raise Exception('Matriz é singular')

def swap_lines(matriz, j, k):
    temp = matriz[j-1]
    matriz[j-1] = matriz[k-1]
    matriz[k-1] = temp

def element(A, i, j):
    return A[i-1][j-1]

def set_element(A, i, j, value):
    A[i-1][j-1] = value

def M_symbols_for_amount(n):
    s = []
    for i in range(n):
        s.append(symbols(f'M{i}'))
    return s


def s_delta(t, h, M, A, B):
    j = 0
    while (t > T[j]):
        j += 1
    i = j - 1

    # SDelta = ..., para t em [ti, ti + 1], sendo Mi, Ai e Bi constantes que dependem de Delta
    resp = (M[i] / 6 * h[i + 1]) * (T[i + 1] - t)**3 + (M[i + 1] / 6 * h[i + 1]) * (t - T[i])**3 + A[i] * (t - T[i]) + B[i]

    return resp


def f(t, h, M, A, B, tau):
    group_number = tau
    distance = 100 - 2 * group_number

    interpolation_interval = binary_search(X, distance)
    print('Intervalo de interpolação:')
    print(interpolation_interval)

    s_t = s_delta(t, h, M, A, B)
    return s_t - distance


def deriv_f(t, h, M, A, B, tau):
    j = 0
    while(t > T[j]):
        j += 1
    i = j - 1

    resp = (-3 * T[i + 1]**2 + 6 * T[i + 1] * t) * (M[i] / 6 * h[i + 1]) + (-6 * t * T[i] + 3 * T[i]**2) * (M[i + 1] / 6 * h[i + 1]) + A[i]
    # ja aplicado, na formula abaixo, a regra da cadeia
    # resp = -3*M[i] / (6*h[i + 1]) * (T[i + 1]**2 + 3*M[i + 1] / (6*h[i + 1]) * (t-T[i])**2 + A[i]

    return resp


def second_deriv_f(t, h, M):
    j = 0
    while(t > T[j]):
        j += 1
    i = j - 1

    resp = (6 * T[i + 1] * M[i] - 6 * T[i] * M[i + 1]) / 6 * h[i + 1]

    return resp


def newton(u_0, epsilon, A, B, M, tau, h):
    u_k = u_0

    while (f(u_k - epsilon, h, M, A, B, tau) * f(u_k + epsilon, h, M, A, B, tau) > 0):
        u_k = u_k - f(u_k, h, M, A, B, tau) / deriv_f(u_k, h, M, A, B, tau)


# SEÇÃO 4 - Algoritmo de busca binária
def binary_search(lista, t):
    m = 0
    M = len(lista) - 1

    while abs(M - m) > 1:
        k = (m + M) // 2
        tk = lista[k]
        if t > tk:
            m = k
        else:
            M = k
    # TODO Deve voltar i = m e não o intervalo
    return [m, M]

# Some function unit tests
def test_swap_lines():
    # Testando duas dimensões
    A = [[0, 1], [2, 3]]
    B = [[2,3],[0,1]]
    swap_lines(A, 1, 2)
    assert(A == B)
    # Testando uma dimensão
    A = [0, 3]
    B = [3, 0]
    swap_lines(A, 2, 1)
    assert(A == B)

def test_find_pivot():
    A = [[0, 1], [2, 3]]
    y = [[2], [3]]
    # print( A )
    # print( y )
    find_pivot(A, y, 1, 2)
    # print( A )
    # print( y )

def test_gauss():
    A = [[2, 3], [1, 1]]
    y = [5, 2]
    answer = gauss(A, y)
    np.testing.assert_array_equal(answer, [1, 1])
    A = [[0, 1], [2, 3]]
    y = [1, 2]
    answer = gauss(A, y)
    np.testing.assert_array_equal(answer, [-0.5, 1])

if __name__ == "__main__":
    # Testes da busca binária
    testlist = [0, 1, 2, 8, 13, 17, 19, 32, 42]
    print(testlist)
    print(binary_search(testlist, 3))
    print(binary_search(testlist, 13))
    # TODO O que deve acontecer se o tempo dado estiver fora dos intervalos?
    print(binary_search(testlist, -1))
    print(binary_search(testlist, 43))
    test_swap_lines()
    test_find_pivot()
    test_gauss()
