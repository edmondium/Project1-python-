import numpy as np
cimport numpy as np
from libc.math cimport exp, M_PI

def ComputeMTDensity(double mu,
                     np.ndarray[np.float64_t, ndim=2] Ek,     # shape (nk, nbands)
                     np.ndarray[np.float64_t, ndim=1] wkp,    # shape (nk,)
                     np.ndarray[np.float64_t, ndim=3] w0,     # shape (nk, nlmax, nbands)
                     np.ndarray[np.float64_t, ndim=3] w1,     # same
                     np.ndarray[np.float64_t, ndim=3] w2,     # same
                     list Psi_l,                              # list of 1D ndarrays
                     list Psip_l,                             # list of 1D ndarrays
                     double beta=50.0):
    """
    Given the coefficients Eqs.70-73 (p.33), computes the valence charge
    with the chemical potential mu. Implements Eq.75 (p.34).
    """

    cdef Py_ssize_t nlmax = w0.shape[1]
    cdef Py_ssize_t nk = Ek.shape[0]
    cdef Py_ssize_t nbands = Ek.shape[1]
    cdef Py_ssize_t l, ik, p, ir
    cdef double x, ferm, dw0, dw1, dw2

    # weights for each l, shape (nlmax, 3)
    cdef np.ndarray[np.float64_t, ndim=2] wgh = np.zeros((nlmax, 3), dtype=np.float64)

    # loop over angular momentum channels
    for l in range(nlmax):
        for ik in range(nk):
            dw0 = dw1 = dw2 = 0.0
            for p in range(nbands):
                x = beta * (Ek[ik, p] - mu)

                if abs(x) < 100.0:
                    ferm = 1.0 / (exp(x) + 1.0)
                else:
                    ferm = 1.0 if x < 0 else 0.0

                dw0 += w0[ik, l, p] * ferm
                dw1 += w1[ik, l, p] * ferm
                dw2 += w2[ik, l, p] * ferm

            wgh[l, 0] += dw0 * wkp[ik]
            wgh[l, 1] += dw1 * wkp[ik]
            wgh[l, 2] += dw2 * wkp[ik]

    # Now compute MTRho
    cdef Py_ssize_t nR = (<np.ndarray> Psi_l[0]).shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] MTRho = np.zeros(nR, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] psi, pspi

    for l in range(nlmax):
        psi = Psi_l[l]
        pspi = Psip_l[l]
        for ir in range(nR):
            MTRho[ir] += (
                wgh[l, 0] * psi[ir] * psi[ir]
                + wgh[l, 1] * psi[ir] * pspi[ir]
                + wgh[l, 2] * pspi[ir] * pspi[ir]
            )

    # multiply by prefactor: 2/(4*pi)
    MTRho *= 2.0 / (4.0 * M_PI)

    return MTRho
