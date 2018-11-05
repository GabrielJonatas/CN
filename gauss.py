import numpy as np

def for_matrix(A):
    n = len(A)
    lastRow = A.shape[0] - 1
    lastCol = A.shape[1] - 1

    for i in range(0, n):
        # Search for maximum in this column
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k

        # Swap maximum row with current row (column by column)
        for k in range(i, n + 1):
            print("\n{}".format(A))
            if k <= lastCol:
                # print('maxRow = {}, max_val = {}, k = {}'.format(maxRow, (A.shape[0] - 1), k))
                tmp = A[maxRow][k]
                A[maxRow][k] = A[i][k]
                A[i][k] = tmp

        # Make all rows below this one 0 in current column
        for k in range(i + 1, n):
            c = -A[k][i] / A[i][i]
            for j in range(i, n + 1):
                # print('k = {}, j = {}'.format(k, j))
                if k <= lastRow and j <= lastCol:
                    if i == j:
                        A[k][j] = 0
                    else:
                        A[k][j] += c * A[i][j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        print("i = {}, n = {}, k = {}".format(i, n, k))
        if i <= lastRow and i <= lastCol and n <= lastCol:
            x[i] = A[i][n] / A[i][i]
            for k in range(i - 1, -1, -1):
                A[k][n] -= A[k][i] * x[i]
    return x


def pprint(A):
    n = len(A)
    for i in range(0, n):
        line = ""
        for j in range(0, n + 1):
            line += str(A[i][j]) + "\t"
            if j == n - 1:
                line += "| "
        print(line)
    print("")


def linearsolver(A, b):
    n = len(A)
    M = A

    i = 0
    for x in M:
        x = np.append(x, b[i])
        i += 1

    for k in range(n):
        for i in range(k, n):
            if abs(M[i][k]) > abs(M[k][k]):
                M[k], M[i] = M[i], M[k]
            else:
                pass

        for j in range(k + 1, n):
            q = float(M[j][k]) / M[k][k]
            for m in range(k, n):
                M[j][m] -= q * M[k][m]

    x = [0 for i in range(n)]

    print('n = ', 3, 'len(M) = ', len(M), 'len(x) = ', len(x))
    x[n - 1] = float(M[n - 1][n]) / M[n - 1][n - 1]
    for i in range(n - 1, -1, -1):
        z = 0
        for j in range(i + 1, n):
            z = z + float(M[i][j]) * x[j]
        x[i] = float(M[i][n] - z) / M[i][i]
    print(x)


if __name__ == "__main__":
    from fractions import Fraction

    n = int(input("\nWhat will be the size of this matrix?\n"))

    A = [[0 for j in range(n + 1)] for i in range(n)]

    # Read input data
    for i in range(0, n):
        line = map(
            Fraction,
            input(
                "\nType the numbers for line {} separated by space:\n".format(i)
            ).split(" "),
        )
        print(list(enumerate(line)))
        for j, el in enumerate(line):
            print("i = {}, j = {}".format(i, j))
            A[i][j] = el

    line = input().split(" ")
    lastLine = list(map(Fraction, line))
    for i in range(0, n):
        A[i][n] = lastLine[i]

    # Print input
    pprint(A)

    # Calculate solution
    x = for_matrix(np.array(A))

    # Print result
    line = "Result:\t"
    for i in range(0, n):
        line += str(x[i]) + "\t"
    print(line)
