# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def cholesky_solve(a,bb):
    """Solve the equation Ax=b where A is a Cholesky-banded matrix.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        :math:`A` in :math:`A x = b`.
    bb : :class:`numpy.ndarray`
        :math:`b` in :math:`A x = b`.

    Returns
    -------
    cholesky_solve : :func:`tuple`
        A tuple containing the status and the result of the solution.  The
        status is always -1.
    """
    import numpy as np
    b = bb.copy()
    bw = a.shape[0]
    n = b.shape[0] - bw
    kn = bw -1
    spot = np.arange(kn,dtype='i4') + 1
    for j in range(n):
        b[j] /= a[0,j]
        b[j+spot] -= b[j]*a[spot,j]
    spot = kn - np.arange(kn,dtype='i4')
    for j in range(n-1,-1,-1):
        b[j] = (b[j] - np.sum(a[spot,j] * b[j+spot]))/a[0,j]
    return (-1,b)
