# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from libc.math cimport NAN

def Numerov(np.ndarray[double, ndim=1] F_in, double dx, double f0=0.0, double f1=1e-3):
    """
    Numerov integration using Cython + typed memoryviews.
    Signature matches your original: Numerov(F, dx, f0=0.0, f1=1e-3)
    Returns a numpy.ndarray[float64] Solution of same length as F_in.
    """
    cdef int Nmax = F_in.shape[0]
    if Nmax <= 0:
        return np.zeros(0, dtype=np.float64)
    if Nmax == 1:
        sol = np.zeros(1, dtype=np.float64)
        sol[0] = f0
        return sol
    # create output array
    cdef np.ndarray[double, ndim=1] Solution = np.zeros(Nmax, dtype=np.float64)

    # create memoryviews for speed
    cdef double[:] F = F_in
    cdef double[:] Sol = Solution

    cdef double h2 = dx * dx
    cdef double h12 = h2 / 12.0

    # initialize
    Sol[0] = f0
    Sol[1] = f1

    cdef double w0 = (1.0 - h12 * F[0]) * Sol[0]
    cdef double Fx = F[1]
    cdef double w1 = (1.0 - h12 * Fx) * Sol[1]
    cdef double Phi = Sol[1]
    cdef double w2
    cdef int i

    for i in range(2, Nmax):
        # w2 = 2*w1 - w0 + h2*Phi*Fx;
        w2 = 2.0 * w1 - w0 + h2 * Phi * Fx
        w0 = w1
        w1 = w2
        Fx = F[i]
        # protect against (1 - h12*Fx) == 0 -> produces inf/nan (keeps parity with original C code)
        if (1.0 - h12 * Fx) == 0.0:
            # choose to set Phi to NAN to indicate problem (consistent with C behaviour of div by zero)
            Phi = NAN
        else:
            Phi = w2 / (1.0 - h12 * Fx)
        Sol[i] = Phi

    return Solution
