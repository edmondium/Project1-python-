# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from libc.math cimport NAN

def NumerovGen(np.ndarray[double, ndim=1] F_in,
               np.ndarray[double, ndim=1] U_in,
               double dx,
               double f0=0.0,
               double f1=1e-3):
    """
    Generalized Numerov integration (inhomogeneous case).
    """
    cdef int Nmax = F_in.shape[0]
    if Nmax <= 0:
        return np.zeros(0, dtype=np.float64)
    if Nmax == 1:
        sol = np.zeros(1, dtype=np.float64)
        sol[0] = f0
        return sol

    cdef np.ndarray[double, ndim=1] Solution = np.zeros(Nmax, dtype=np.float64)

    # memoryviews
    cdef double[:] F = F_in
    cdef double[:] U = U_in
    cdef double[:] Sol = Solution

    cdef double h2 = dx * dx
    cdef double h12 = h2 / 12.0

    # initial conditions
    Sol[0] = f0
    Sol[1] = f1

    cdef double w0 = Sol[0] * (1.0 - h12 * F[0]) - h12 * U[0]
    cdef double Fx = F[1]
    cdef double Ux = U[1]
    cdef double w1 = Sol[1] * (1.0 - h12 * Fx) - h12 * Ux
    cdef double Phi = Sol[1]

    cdef double w2
    cdef int i

    for i in range(2, Nmax):
        w2 = 2.0 * w1 - w0 + h2 * (Phi * Fx + Ux)
        w0 = w1
        w1 = w2
        Fx = F[i]
        Ux = U[i]
        if (1.0 - h12 * Fx) == 0.0:
            Phi = NAN  # division by zero guard
        else:
            Phi = (w2 + h12 * Ux) / (1.0 - h12 * Fx)
        Sol[i] = Phi

    return Solution
