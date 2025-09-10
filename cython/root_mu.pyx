import numpy as np
cimport numpy as np
from libc.math cimport exp

def rootChemicalPotential(double mu,
                          np.ndarray[np.float64_t, ndim=2] Ek,
                          np.ndarray[np.float64_t, ndim=1] wkp,
                          double Zval,
                          double beta=50.0):
    """
    Computes valence density to find root for the chemical potential
    """

    cdef Py_ssize_t ik, p
    cdef Py_ssize_t nk = Ek.shape[0]
    cdef Py_ssize_t npbands = Ek.shape[1]
    cdef double Zt = 0.0
    cdef double x, ferm

    for ik in range(nk):
        for p in range(npbands):
            x = beta * (Ek[ik, p] - mu)

            if abs(x) < 100.0:
                ferm = 1.0 / (exp(x) + 1.0)
            else:
                ferm = 1.0 if x < 0 else 0.0

            Zt += wkp[ik] * ferm

    return 2.0 * Zt - Zval
