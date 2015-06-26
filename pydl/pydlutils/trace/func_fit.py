# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def func_fit(x,y,ncoeff,invvar=None,function_name='legendre',ia=None,inputans=None,inputfunc=None):
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
    inputfunc : array-like, optional
        Multiply the function fit by these values.

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
    ngood = len(igood)
    res = np.zeros((ncoeff,),dtype=x.dtype)
    yfit = np.zeros(x.shape,dtype=x.dtype)
    if ngood == 0:
        pass
    elif ngood == 1:
        res[0] = y[igood[0]]
        yfit += y[igood[0]]
    else:
        ncfit = min(ngood,ncoeff)
        function_map = {
            'legendre':flegendre,
            'flegendre':flegendre,
            'chebyshev':fchebyshev,
            'fchebyshev':fchebyshev,
            'chebyshev_split':fchebyshev_split,
            'fchebyshev_split':fchebyshev_split,
            'poly':fpoly,
            'fpoly':fpoly
            }
        try:
            legarr = function_map[function_name](x,ncfit)
        except KeyError:
            raise KeyError('Unknown function type: {0}'.format(function_name))
        if inputfunc is not None:
            if inputfunc.shape != x.shape:
                raise ValueError('Dimensions of X and inputfunc do not agree!')
            legarr *= np.tile(inputfunc,ncfit).reshape(ncfit,x.shape[0])
        yfix = np.zeros(x.shape,dtype=x.dtype)
        nonfix = ia[0:ncfit].nonzero()[0]
        nparams = len(nonfix)
        fixed = (~ia[0:ncfit]).nonzero()[0]
        if len(fixed) > 0:
            yfix = np.dot(legarr.T,inputans * (1 - ia))
            ysub = y - yfix
            finalarr = legarr[nonfix,:]
        else:
            finalarr = legarr
            ysub = y
        # extra2 = finalarr * np.outer(np.ones((nparams,),dtype=x.dtype),(invvar > 0))
        extra2 = finalarr * np.outer(np.ones((nparams,),dtype=x.dtype),invvar)
        alpha = np.dot(finalarr,extra2.T)
        # assert alpha.dtype == x.dtype
        if nparams > 1:
            # beta = np.dot(ysub * (invvar > 0), finalarr.T)
            beta = np.dot(ysub * invvar, finalarr.T)
            assert beta.dtype == x.dtype
            # uu,ww,vv = np.linalg.svd(alpha,full_matrices=False)
            res[nonfix] = np.linalg.solve(alpha,beta)
        else:
            # res[nonfix] = (ysub * (invvar > 0) * finalarr).sum()/alpha
            res[nonfix] = (ysub * invvar * finalarr).sum()/alpha
        if len(fixed) > 0:
            res[fixed] = inputans[fixed]
        yfit = np.dot(legarr.T,res[0:ncfit])
    return (res,yfit)
