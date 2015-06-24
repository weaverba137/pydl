# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def func_fit(x,y,ncoeff,invvar=None,function_name='legendre',ia=None,inputans=None):
    """Fit X, Y positions to a functional form.

    Parameters
    ----------
    x : array-like
        X values (independent variable).
    y : array-like
        Y values (dependent variable).
    ncoeff : int
        Number of coefficients to fit.
    invvar : array-like, optional
        Weight values; inverse variance.
    function_name : str, optional
        Function name, default 'legendre'.
    ia : array-like, optional
        An array of bool of length `ncoeff` specifying free (``True``) and
        fixed (``False``) parameters.
    inputans : array-like, optional
        An array of values of length `ncoeff` specifying the values of
        the fixed parameters.

    Returns
    -------
    func_fit : tuple of array-like
        Fit coefficients, length `ncoeff`; fitted values.
    """
    import numpy as np
    from . import fchebyshev, fchebyshev_split, fpoly
    from ...goddard.math import flegendre
    if x.shape != y.shape:
        raise ValueError('Dimensions of X and Y do not agree!')
    if invvar is None:
        invvar = np.ones(x.shape,dtype=x.dtype)
    else:
        if invvar.shape != x.shape:
            raise ValueError('Dimensions of X and invvar do not agree!')
    if ia is None:
        ia = np.ones((ncoeff,),dtype=np.bool)
    if not ia.all():
        if inputans is None:
            inputans = np.zeros((ncoeff,),dtype=x.dtype)
    #
    # Select unmasked points
    #
    igood = (invvar > 0).nonzero()[0]
    ngood = igood.sum()
    res = np.zeros((ncoeff,),dtype=x.dtype)
    yfit = np.zeros(x.shape,dtype=x.dtype)
    if ngood == 0:
        pass
    elif ngood == 1:
        res[0] = y[igood[0]]
        yfit += y[igood[0]]
    else:
        pass
    return (res,yfit)
