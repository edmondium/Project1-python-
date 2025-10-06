# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np

# Tell Cython the expected array types
def CRHS(double E, int l,
         np.ndarray[double, ndim=1] R,
         np.ndarray[double, ndim=1] Veff):
    """
    RHS for solving the Schroedinger equations by Numerov.
    Cython implementation replacing scipy.weave.
    """
    cdef int N = R.shape[0]
    cdef np.ndarray[double, ndim=1] RHS = np.zeros(N, dtype=np.float64)

    cdef int i
    cdef double Ri
    for i in range(N):
        Ri = R[i]
        RHS[i] = 2.0 * (-E + 0.5 * l * (l + 1) / (Ri * Ri) + Veff[i])

    return RHS
