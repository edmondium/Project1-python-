# cython: boundscheck=False, wraparound=False, cdivision=True, wraparound=False
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

def compute_argums(np.ndarray[double, ndim=2] Km_in,
                   np.ndarray[double, ndim=1] k_in):
    """
    Compute qv, qvs, argums from reciprocal vectors.
    - Km_in: shape (N,3)
    - k_in: shape (3,)
    Returns (qv, qvs, argums)
    """
    cdef int N = Km_in.shape[0]
    cdef np.ndarray[double, ndim=2] qv = np.zeros((N, 3), dtype=np.float64)
    cdef np.ndarray[double, ndim=1] qvs = np.zeros(N, dtype=np.float64)
    cdef np.ndarray[double, ndim=2] argums = np.zeros((N, N), dtype=np.float64)

    # memoryviews for speed
    cdef double[:, :] Km = Km_in
    cdef double[:] k = k_in
    cdef double[:, :] qv_mv = qv
    cdef double[:] qvs_mv = qvs
    cdef double[:, :] argums_mv = argums

    cdef int iK, jK, i
    cdef double qvqv

    # compute qv and qvs
    for iK in range(N):
        for i in range(3):
            qv_mv[iK, i] = Km[iK, i] + k[i]
        qvs_mv[iK] = sqrt(qv_mv[iK, 0]*qv_mv[iK, 0] +
                          qv_mv[iK, 1]*qv_mv[iK, 1] +
                          qv_mv[iK, 2]*qv_mv[iK, 2])

    # compute argums
    for iK in range(N):
        for jK in range(N):
            qvqv = 0.0
            for i in range(3):
                qvqv += qv_mv[iK, i] * qv_mv[jK, i]
            if qvs_mv[iK] * qvs_mv[jK] == 0.0:
                argums_mv[iK, jK] = 1.0
            else:
                argums_mv[iK, jK] = qvqv / (qvs_mv[iK] * qvs_mv[jK])

    return qv, qvs, argums
