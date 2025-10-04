# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np

def compute_HO(np.ndarray[double, ndim=2] Olap_I_in,   # (nK, nK)
               np.ndarray[double, ndim=1] qvs_in,     # (nK,)
               np.ndarray[double, ndim=3] Leg_in,     # (nK, nK, nL)
               np.ndarray[double, ndim=1] PP_in,      # (nL,)
               np.ndarray[double, ndim=2] a_lk_in,    # (nK, nL)  (note: indexing a_lk[iK, l])
               np.ndarray[double, ndim=2] b_lk_in,    # (nK, nL)
               np.ndarray[double, ndim=1] Enu_in,     # (nL,)
               double VKSi,
               double RMuffinTin,
               double Vol):
    """
    Compute:
      Olap (nK,nK), Ham (nK,nK),
      WK0, WK1, WK2 (each (nL, nK, nK))
    Based on your original C++/weave code.
    """

    cdef int nK = Olap_I_in.shape[0]
    cdef int nK2 = Olap_I_in.shape[1]
    if nK != nK2:
        raise ValueError("Olap_I must be square (nK,nK)")
    cdef int nL = Enu_in.shape[0]

    # allocate outputs
    cdef np.ndarray[double, ndim=2] Olap = np.zeros((nK, nK), dtype=np.float64)
    cdef np.ndarray[double, ndim=2] Ham  = np.zeros((nK, nK), dtype=np.float64)
    cdef np.ndarray[double, ndim=3] WK0  = np.zeros((nL, nK, nK), dtype=np.float64)
    cdef np.ndarray[double, ndim=3] WK1  = np.zeros((nL, nK, nK), dtype=np.float64)
    cdef np.ndarray[double, ndim=3] WK2  = np.zeros((nL, nK, nK), dtype=np.float64)

    # memoryviews for inputs
    cdef double[:, :] Olap_I = Olap_I_in
    cdef double[:] qvs = qvs_in
    cdef double[:, :, :] Leg = Leg_in
    cdef double[:] PP = PP_in
    cdef double[:, :] a_lk = a_lk_in
    cdef double[:, :] b_lk = b_lk_in
    cdef double[:] Enu = Enu_in

    # memoryviews for outputs
    cdef double[:, :] Olap_mv = Olap
    cdef double[:, :] Ham_mv = Ham
    cdef double[:, :, :] WK0_mv = WK0
    cdef double[:, :, :] WK1_mv = WK1
    cdef double[:, :, :] WK2_mv = WK2

    # local variables
    cdef int iK, jK, l
    cdef double olapMT, hamMT, Plk
    cdef double a_a, a_b, b_a, b_b, olap_local
    cdef double factor = np.pi * RMuffinTin * RMuffinTin * RMuffinTin * RMuffinTin / Vol
    # C0n(l) = (2*l + 1) * pi * RMuffinTin**4 / Vol

    for iK in range(nK):
        for jK in range(nK):
            olapMT = 0.0
            hamMT = 0.0
            for l in range(nL):
                Plk = (2.0 * l + 1.0) * factor * Leg[iK, jK, l]
                # note: a_lk and b_lk indexing matches original a_lk(iK,l)
                a_a = a_lk[iK, l] * a_lk[jK, l]
                a_b = a_lk[iK, l] * b_lk[jK, l]
                b_a = b_lk[iK, l] * a_lk[jK, l]
                b_b = b_lk[iK, l] * b_lk[jK, l]

                olap_local = a_a + b_b * PP[l]
                olapMT += Plk * olap_local
                hamMT  += Plk * (0.5 * (b_a + a_b) + olap_local * Enu[l])

                WK0_mv[l, iK, jK] = Plk * a_a
                WK1_mv[l, iK, jK] = Plk * (b_a + a_b)
                WK2_mv[l, iK, jK] = Plk * b_b

            Olap_mv[iK, jK] = olapMT + Olap_I[iK, jK]
            Ham_mv[iK, jK] = (0.25 * (qvs[iK] * qvs[iK] + qvs[jK] * qvs[jK]) + VKSi) * Olap_I[iK, jK] + hamMT

    return Olap, Ham, WK0, WK1, WK2
