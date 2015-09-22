# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def cholesky_band(l,mininf=0.0,verbose=False):
    """Compute Cholesky decomposition of banded matrix.

    Parameters
    ----------
    l : :class:`numpy.ndarray`
        A matrix on which to perform the Cholesky decomposition.
    mininf : :class:`float`, optional
        Entries in the `l` matrix are considered negative if they are less
        than this value (default 0.0).
    verbose : :class:`bool`, optional
        If set to ``True``, print some debugging information.

    Returns
    -------
    cholesky_band : :func:`tuple`
        If problems were detected, the first item will be the index or
        indexes where the problem was detected, and the second item will simply
        be the input matrix.  If no problems were detected, the first item
        will be -1, and the second item will be the Cholesky decomposition.
    """
    import numpy as np
    lower = l.copy()
    bw,nn = lower.shape
    n = nn - bw
    if verbose:
        print(lower[0,0:n])
    negative = lower[0,0:n] <= mininf
    if negative.any() or not np.all(np.isfinite(lower)):
        if verbose:
            print('Bad entries.')
            print(negative.nonzero()[0])
        return (negative.nonzero()[0],l)
    kn = bw - 1
    spot = np.arange(kn,dtype='i4') + 1
    bi = np.arange(kn,dtype='i4')
    for i in range(1,kn):
        bi = np.append(bi, np.arange(kn-i,dtype='i4') + (kn+1)*i)
    for j in range(n):
        lower[0,j] = np.sqrt(lower[0,j])
        lower[spot,j] /= lower[0,j]
        x = lower[spot,j]
        if not np.all(np.isfinite(x)):
            if verbose:
                print('NaN found in cholesky_band.')
            return (j,l)
        hmm = np.outer(x,x)
        here = bi+(j+1)*bw
        lower.T.flat[here] -= hmm.flat[bi]
    return (-1,lower)
