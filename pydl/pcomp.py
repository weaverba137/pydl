# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def pcomp(x,**kwargs):
    """Replicates the IDL PCOMP() function.

    Parameters
    ----------
    x : array_like
        A 2-D arrray.

    Returns
    -------
    pcomp : dict

    Notes
    -----

    References
    ----------
    http://www.exelisvis.com/docs/PCOMP.html

    Examples
    --------
    """
    from numpy import corrcoef, cov, dot, sqrt, tile
    from scipy.linalg import eigh
    if x.ndim != 2:
        raise ValueError('Input array must be two-dimensional')
    no,nv = x.shape
    if 'standardize' in kwargs:
        xstd = x - tile(x.mean(0),no).reshape(x.shape)
        s = tile(xstd.std(0),no).reshape(x.shape)
        array = xstd/s
    else:
        array = x
    if 'covariance' in kwargs:
        c = cov(array,rowvar=0)
    else:
        c = corrcoef(array,rowvar=0)
    #
    # eigh is used for symmetric matrices
    #
    evals, evecs = eigh(c)
    #
    # Sort eigenvalues in descending order
    #
    ie = evals.argsort()[::-1]
    evals = evals[ie]
    evecs = evecs[:,ie]
    #
    # If necessary, add code to fix the signs of the eigenvectors.
    # http://www3.interscience.wiley.com/journal/117912150/abstract
    #
    normevecs = evecs * tile(sqrt(evals),nv).reshape(nv,nv)
    variances = evals/c.trace()
    derived_data = dot(array,normevecs)
    if 'standardize' in kwargs:
        derived_data += xstd
    return {'derived':derived_data,'coefficients':normevecs,
        'variance':variances,'eigenvalues':evals}

