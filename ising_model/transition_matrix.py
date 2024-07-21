import numpy as np

np.set_printoptions(precision=7)

def self_matrix_product(M, n):
    for i in range(n-1):
        M = M @ M
    return M


def self_matrix_product2(M, n):
    if n - 1 == 0:
        return M
    return M @ self_matrix_product2(M, n-1)


A = np.array([[1/2, 1/2], [3/4, 1/4]])
B = self_matrix_product2(A, 3)
C = self_matrix_product2(A, 9)
D = self_matrix_product2(A, 10)
print(A)
print(B)
print(C)
print(D)