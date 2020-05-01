import math
import numpy as np


def jacobi(A, d, mi, lambda_array, n_max=100):
    # default
    tolerance = 10**-8
    epsilon = tolerance + 1
    k = 0
    x_0 = [0 for i in range(29)]
    x_old = list(np.copy(x_0))
    x_new = list(np.copy(x_0))
    n = len(x_0)

    sigma = np.zeros(len(A))
    for i in range(len(A)):
        # print(f'i = {i}\nsigma[i] = {sigma[i]}\nA[i][i] = {A[i][i]}')
        sigma[i] = (1 / A[i][i]) * (sum(A[i]) - A[i][i])
    sigma = max(sigma)
    # print(f'sigma = {sigma}')

    while epsilon > tolerance and k < n_max:
        # print(f'[JACOBI] Iteration {k}'
        if k == 1:

            norm = infinity_norm(x_new, x_0)
            print(f'Norma infinita: {norm}'.format(norm))

            minimum_n_decimal = math.log((tolerance * (1 - sigma)) / (norm)) / math.log(sigma)
            print(f'Minimo numero de iterações, decimal: {minimum_n_decimal}')
            n_max = math.ceil(minimum_n_decimal)
            print(f'Minimo numero de iterações, arredondado: {n_max}')

        x_new[0] = .5 * (d[0] - lambda_array[0] * x_old[1])

        for i in range(1, n - 1):
            x_new[i] = .5 * (d[i] - mi[i] * x_old[i - 1] - lambda_array[i] * x_old[i + 1])

        x_new[n - 1] = .5 * (d[len(d) - 1] - mi[len(mi) - 1] * x_old[n - 1])

        epsilon = infinity_norm(x_new, x_old)

        x_old = list(np.copy(x_new))
        k += 1

    return x_new


def subtract(x1, x2):
    return x1 - x2


def infinity_norm(v1, v2):
    new_minus_old = list(map(subtract, v1, v2))
    new_minus_old = [-x if x < 0 else x for x in new_minus_old]
    return max(new_minus_old)
