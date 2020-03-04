import numpy as np

# def triag(n, arr1, arr2, free_coeffs_A):
#     A = np.zeros((2*n, 2*n)) - np.diag(np.ones(2*n-1), k=-1) +\
#         np.diag(np.ones(2*n-1), k=1)
#     A = np.append(np.delete(A, 0, axis=0), [arr2], axis=0)
#     A[-2, :] = arr1
#     b = np.zeros(2*n)
#     for i in range(0, len(A)-2, 2):
#         b[i] = np.real(free_coeffs_A[i//2])
#         b[i+1] = np.imag(free_coeffs_A[i//2])
#     return np.linalg.solve(A, b)

def toDot(r, F):
    return np.sum(np.real(r)*np.real(F) + np.imag(r)*np.imag(F))

def toCross(r, F):
    return np.sum(-np.real(F)*np.imag(r) + np.real(r)*np.imag(F))

def moment(rs):
    return np.array([-np.imag(rs), np.real(rs)])

def mat(a_List, newshape):
    A = [nested_element for element in a_List for nested_element in element]
    return np.reshape(A, newshape)