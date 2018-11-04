import numpy as np

T = np.arange(0, 1.1, .1, dtype=float)
X = [2.00, 1.81, 1.64, 1.49, 1.36, 1.25, 1.16, 1.09, 1.04, 1.01, 1]

def monta_h(v):
	n = len(v)
	# Instancia h como vetor vazio com tamanho n - 1
	h = np.empty(n - 1)
	for i in range(n - 1):
		h[i] = v[i+1] - v[i]
	return h

def monta_matriz(n):
	# Instancia matrix nxn
	# preenchida com 0
	A = np.full((n, n), 0.0, dtype=float)

	for i in range(n):
		A[i][i] = 2
		A[i][i+1] = .5
		A[i+1][i] = .5
	
	return A

# x00 = 2a DERIVADA DE x0
# xnn = 2a DERIVADA DE xn
# calcule x00 & xnn A PARTIR DO DEFINIDO NO SLIDE, TEM FORMULA DE APROXIMACAO PRA AMBOS
def monta_vetor(n, x00, xnn):
	# Instancia matrix nxn
	# preenchida com 0
	d = np.full((n, n), 0.0, dtype=float)

	d[1,1] = 2*x00
	d[n+1, 1] = 2*xnn
	
	for i in range(2, n+1):
		d[i, 1] = 6*(2*(X[i+1] - X[i]) - 2*X*(X[i] - X[i-1]))
	
	return d

# VETOR DE VALORES M 
def M(n, x00, xnn):
	A = monta_matriz(n)
	d = monta_vetor(n, x00, xnn)
	# implementar Gauss e chamar a funcao do dito cujo
	return gauss(A, d)

# solucao de sistema por Gauss
def gauss(A, d):
	# DO THISSSSSSSSSSSSSS
	# https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
	return

def s_delta(t, h, M, A, B):
	j = 0
	while (t > T[j]):
		j += 1
	i = j - 1

	# SDelta = ..., para t em [ti, ti+1], sendo Mi, Ai e Bi constantes que dependem de Delta
	resp = (M[i] / 6*h[i+1])*(T[i+1] - t)^3 + (M[i+1] / 6*h[i+1])*(t - T[i])^3 + A[i]*(t - T[i]) + B[i]
	
	return resp

# def deriv_f(t, h, M, A, B, tau):
# 	j = 0
# 	while(t > T[j]):
# 		j += 1
# 	i = j - 1
	
# 	# ja aplicado, na formula abaixo, a regra da cadeia
# 	resp = -3*M[i] / (6*h[i+1]) * (T[i+1]^2 + 3*M[i+1] / (6*h[i+1]) * (t-T[i])^2 + A[i]
	

def newton(u_0, f, DerivF, epsilon, A, B, tau, h):
	u_k = u_0
	
	# passe as variaveis ", h, M, A, B, tau" à f() abaixo OU defina-as como variáveis globais/acessíveis para f()
	while (f(u_k - epsilon, h, M, A, B, tau) * f(u_k + epsilon, h, M, A, B, tau) > 0):
		u_k = u_k - f(u_k)/DerivF(u_k)

# SEÇÃO 4 - Algoritmo de busca binária
def buscaBinaria(lista, t):
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

# Testes da busca binária
testlist = [0, 1, 2, 8, 13, 17, 19, 32, 42,]
print(testlist)
print(buscaBinaria(testlist, 3))
print(buscaBinaria(testlist, 13))
# print("================================")
# print(xbuscaBinaria(testlist, 3))
# print(xbuscaBinaria(testlist, 13))