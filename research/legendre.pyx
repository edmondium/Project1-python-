# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np

def compute_legendre(np.ndarray[double, ndim=2] argums_in,
                     int lmax):
    """
    Compute Legendre polynomials P_l(x) for all entries in argums.
    - argums_in: shape (N, N), containing cos(theta) values
    - lmax: maximum order
    Returns: Leg (N, N, lmax+1)
    """

    cdef int N1 = argums_in.shape[0]
    cdef int N2 = argums_in.shape[1]

    cdef np.ndarray[double, ndim=3] Leg = np.zeros((N1, N2, lmax+1),
                                                   dtype=np.float64)

    # memoryviews
    cdef double[:, :] argums = argums_in
    cdef double[:, :, :] Leg_mv = Leg

    cdef int iK, jK, l, i
    cdef double x, x2
    cdef double p0, p1, p2

    for iK in range(N1):
        for jK in range(N2):
            x = argums[iK, jK]
            x2 = x * x

            Leg_mv[iK, jK, 0] = 1.0
            if lmax >= 1:
                Leg_mv[iK, jK, 1] = x
            if lmax >= 2:
                Leg_mv[iK, jK, 2] = 1.5 * x2 - 0.5
            if lmax >= 3:
                Leg_mv[iK, jK, 3] = x * (2.5 * x2 - 1.5)
            if lmax >= 4:
                Leg_mv[iK, jK, 4] = 0.375 * (1 - 10*x2*(1 - 1.1666666666666667*x2))
            if lmax >= 5:
                Leg_mv[iK, jK, 5] = 1.875 * x * (1 - 4.66666666666666667*x2*(1 - 0.9*x2))

            if lmax >= 6:
                # initialize recurrence
                p0 = 0.375 * (1 - 10*x2*(1 - 1.1666666666666667*x2))
                p1 = 1.875 * x * (1 - 4.66666666666666667*x2*(1 - 0.9*x2))
                p2 = 0.0
                for l in range(6, lmax+1):
                    p2 = ((2*l - 1)*x*p1 - (l - 1)*p0) / l
                    p0 = p1
                    p1 = p2
                Leg_mv[iK, jK, lmax] = p2

    return Leg
