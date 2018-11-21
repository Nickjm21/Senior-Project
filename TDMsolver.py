import numpy as np


# Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_
    (Thomas_algorithm)
    '''
    nf = len(d)  # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))  # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]

    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc


A = np.array([[1, 4-2j, 0, 0], [2+1j, 3+1j, 12+2j, 0],
              [0, 7+3j, 8-1j, 4+1j], [0, 0, 10-8j, -4j]], dtype=complex)

a = np.array([2.+1.j, 7.+3.j, 10.-8.j], dtype=complex)
b = np.array([1., 3.+1.j, 8.-1.j, -4.j], dtype=complex)
c = np.array([4.-2.j, 12.+2.j, 4.+1.j], dtype=complex)
d = np.array([7.-1.j, 18.+4.j, -7.5j, 4.-7.j], dtype=complex)

# TDMAsolver(a, b, c, d)
# compare against numpy linear algebra library
# np.linalg.solve(A, d)
