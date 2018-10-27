def xbuscaBinaria(lista, item):
	m = 0
	M = len(lista)-1
	found = False

	while m<=M and not found:
		k = (m + M)//2
		print("k = {} para m = {} & M = {}".format(k, m, M))
		if lista[k] == item:
			found = True
		else:
			if item < lista[k]:
				M = k-1
			else:
				m = k+1
	print([m, M])
	return found


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