import numpy as np
from math import log

def jacobi(A, d, mi, lambda_array):
        # default
        n_max = 100
        tolerance = np.e**-4
        epsilon = 10**-8
        k = 0
        x_0 = [0 for i in range(0, 30)]
        x_old = x_0
        x_new = x_0
        n = len(x_0)

        sigma = np.zeros(len(A))
        for i in range(len(A)):
            # print('i = {}\nsigma[i] = {}\nA[i][i] = {}'.format(i, sigma[i], A[i][i]))
            sigma[i] = (1/A[i][i]) * (sum(A[i]) - A[i][i])


        while epsilon > tolerance and k < n_max:
            print('[JACOBI] Iteration {}'.format(k))

            x_new[0] = .5 * (d[0] - lambda_array[0]*x_old[1])

            for i in range(1, n-1):
                x_new[i] = .5 * (d[i] - mi[i]*x_old[i-1] - lambda_array[i]*x_old[i+1])
            
            x_new[n] = .5 * (d[n] - mi[n]*x_old[n-1])


            subtract = lambda x1, x2: x1 - x2
            new_minus_old = list(map(subtract, x_new, x_old))
            new_minus_old = [ -x if x < 0 else x for x in new_minus_old ]
            epsilon = max(new_minus_old)

            x_old = x_new
            k += 1
        
        return x_new