import numpy as np
from sympy import Matrix, linsolve, solve_linear_system, symbols

T = np.arange(0, 30)
epsilon = 10**-3
X = [
	0, 1, 2.4, 4.1, 6, 8.2, 10.6, 13.4, 16.4, 19.7,
	23.3, 27, 31.2, 35.5, 40.1, 45, 50.2, 55.6, 61.3, 67.3,
	73.6, 80.1, 86.9, 94, 101.3, 109, 116.9, 125, 133.4, 142.1
]

x00 = lambda T=T: ( (X[2] - X[1])/(T[2] - T[1]) - (X[1] - X[0])/(T[1] - T[0]) )/(T[2] - T[0])
xnn = lambda T, n: ( (X[n] - X[n-1])/(T[n] - T[n-1]) - (X[n-1] - X[n-2])/(T[n-1] - T[n-2]) )/(T[n] - T[n-2])

def h(T):
	n = len(T)
	h = np.zeros(n - 1)

	# Valores de h são as diferenças entre valores de t
	for i in range(n - 1):
		h[i] = T[i+1] - T[i]
	return h

def mi_array(h):
	n = len(h) - 1
	mi = np.zeros(n)

	for i in range(1, n):
		mi[i-1] = h[i]/(h[i]+h[i+1])

	return mi

def lambda_array(h):
	n = len(h) - 1
	lambda_array = np.zeros(n)

	for i in range(n):
		lambda_array[i] = 0 if i == 0 else h[i+1]/(h[i]+h[i+1])

	return lambda_array

def make_matrix(n, mi_array, lambda_array):
	# Instancia matrix nxn
	A = np.zeros((n, n))

	for i in range(n):
		A[i][i] = 2
		if i != (n - 1):
			A[i][i+1] = lambda_array[i]
			A[i+1][i] = mi_array[i]
	
	return A

def make_d(h, x00, xnn):
	n = len(h)
	d = np.zeros(n)
	
	d[0] = 2*x00
	d[n-1] = 2*xnn

	for i in range(1, n-1):
		d[i] = (6/(h[i] + h[i+1]))*( (X[i+1]-X[i]/h[i+1]) - (X[i] - X[i-1])/h[i] )

	return d

def B(M, h, n):
	B = np.zeros(n)

	for i in range(n - 1):
		B[i] = X[i] - (M[i]/6) * h[i+1]**2

	return B

def A(M, h, n):
	A = np.full(n, 0.0, dtype=float)

	for i in range(n - 1):
		A[i] = (X[i+1] - X[i])/h[i+1] - ((M[i+1] - M[i])/6)*h[i+1]
		print('A[{}] = {}'.format(i, A[i]))

	return A


# VETOR DE VALORES M 
def M_gauss(n, x00, xnn, mi, lambda_array, h):
	A = make_matrix(n, mi, lambda_array)
	d = make_d(h, x00, xnn)
	return gauss(A, d)

def gauss(A, d):
	print('A =\n{}\nd =\n{}'.format(A, d))

	C = np.zeros((len(A), len(A[0]) + 1))
	print('len(A) = {}, len(A[0]) = {}'.format(len(A), len(A[0])))
	print('len(C) = {}, len(C[0]) = {}'.format(len(C), len(C[0])))
	print('len(d) = {}, d = {}'.format(len(d), d))

	# Preenche matriz C com o sistema completo: Matriz A + Vetor d	
	for i in range(len(A)):
		for j in range(len(A[0]) + 1):
			C[i][j] = A[i][j] if j < len(A[0]) else d[i]
	system = Matrix(C)

	symbols = M_symbols_for_amount(len(A))
	return solve_linear_system(system, *symbols)

def M_symbols_for_amount(n):
	s = []
	for i in range(n):
		s.append(symbols('M{}'.format(i)))
	return s

def s_delta(t, h, M, A, B):
	j = 0
	while (t > T[j]):
		j += 1
	i = j - 1

	# SDelta = ..., para t em [ti, ti+1], sendo Mi, Ai e Bi constantes que dependem de Delta
	resp = (M[i] / 6*h[i+1])*(T[i+1] - t)**3 + (M[i+1] / 6*h[i+1])*(t - T[i])**3 + A[i]*(t - T[i]) + B[i]
	
	return resp

# def deriv_f(t, h, M, A, B, tau):
# 	j = 0
# 	while(t > T[j]):
# 		j += 1
# 	i = j - 1
	
# 	# ja aplicado, na formula abaixo, a regra da cadeia
# 	resp = -3*M[i] / (6*h[i+1]) * (T[i+1]**2 + 3*M[i+1] / (6*h[i+1]) * (t-T[i])**2 + A[i]
	

def newton(u_0, f, DerivF, epsilon, A, B, M, tau, h):
	u_k = u_0
	
	# passe as variaveis "h, M, A, B, tau" à f() abaixo OU defina-as como variáveis globais/acessíveis para f()
	while (f(u_k - epsilon, h, M, A, B, tau) * f(u_k + epsilon, h, M, A, B, tau) > 0):
		u_k = u_k - f(u_k)/DerivF(u_k)

# SEÇÃO 4 - Algoritmo de busca binária
def binary_search(lista, t):
	print("t = {}".format(t))
	m = 0
	M = len(lista)-1

	while abs(M - m) > 1:
		k = (m + M)//2
		tk = lista[k]
		print("k = {} para m = {} & M = {} [tk = {}]".format(k, m, M, tk))
		if t > tk:
			m = k
		if t <= tk:
			M = k
	return [m, M]

if __name__ == "__main__":
	# Testes da busca binária
	testlist = [0, 1, 2, 8, 13, 17, 19, 32, 42,]
	print(testlist)
	print(binary_search(testlist, 3))
	print(binary_search(testlist, 13))