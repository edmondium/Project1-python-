# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np

def SolvePoisson(double Zq,
                 np.ndarray[double, ndim=1] R_in,
                 np.ndarray[double, ndim=1] rho_in):
    """
    Given the input density rho, calculates the Hartree potential.
    Boundary conditions: U(0)=0 and U(S)=Zq (shifted later).
    """

    cdef int Nmax = R_in.shape[0]
    cdef np.ndarray[double, ndim=1] Solution = np.zeros(Nmax, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] U = np.zeros(Nmax, dtype=np.float64)

    # typed memoryviews
    cdef double[:] R = R_in
    cdef double[:] rho = rho_in
    cdef double[:] Sol = Solution
    cdef double[:] Uv = U

    cdef int i
    cdef double r, dx, h2, h12
    cdef double w0, w1, w2, Phi, Ux

    # Build U = -4π r * rho[i]
    for i in range(Nmax):
        r = R[i]
        Uv[i] = -4.0 * np.pi * r * rho[i]

    dx = (R[Nmax - 1] - R[0]) / (Nmax - 1.0)
    h2 = dx * dx
    h12 = h2 / 12.0

    # Boundary conditions
    Sol[0] = 0.0
    Sol[1] = (R[1] - R[0])  # Boundary condition for U_H = V_H / r

    # Initialize Numerov inhomogeneous recurrence
    w0 = Sol[0] - h12 * Uv[0]
    Ux = Uv[1]
    w1 = Sol[1] - h12 * Ux
    Phi = Sol[1]

    for i in range(2, Nmax):
        w2 = 2.0 * w1 - w0 + h2 * Ux
        w0 = w1
        w1 = w2
        Ux = Uv[i]
        Phi = w2 + h12 * Ux
        Sol[i] = Phi

    # Add homogeneous solution to satisfy boundary condition at infinity
    cdef double alpha = (Zq - Sol[Nmax - 1]) / R[Nmax - 1]
    for i in range(Nmax):
        Sol[i] += alpha * R[i]

    return Solution
