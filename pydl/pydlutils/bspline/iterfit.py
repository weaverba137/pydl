# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def iterfit(xdata,ydata,invvar=None,upper=5,lower=5,x2=None,maxiter=10,**kwargs):
    """Iteratively fit a b-spline set to data, with rejection.

    Parameters
    ----------
    xdata : :class:`numpy.ndarray`
        Independent variable.
    ydata : :class:`numpy.ndarray`
        Dependent variable.
    invvar : :class:`numpy.ndarray`
        Inverse variance of `ydata`.  If not set, it will be calculated based
        on the standard deviation.
    upper : :class:`int` or :class:`float`
        Upper rejection threshold in units of sigma, defaults to 5 sigma.
    lower : :class:`int` or :class:`float`
        Lower rejection threshold in units of sigma, defaults to 5 sigma.
    x2 : :class:`numpy.ndarray`, optional
        Orthogonal dependent variable for 2d fits.
    maxiter : :class:`int`, optional
        Maximum number of rejection iterations, default 10.  Set this to
        zero to disable rejection.

    Returns
    -------
    iterfit : :func:`tuple`
        A tuple containing the fitted bspline object and an output mask.
    """
    import numpy as np
    from .bspline import bspline
    from ..math import djs_reject
    nx = xdata.size
    if ydata.size != nx:
        raise ValueError('Dimensions of xdata and ydata do not agree.')
    if invvar is not None:
        if invvar.size != nx:
            raise ValueError('Dimensions of xdata and invvar do not agree.')
    else:
        #
        # This correction to the variance makes it the same
        # as IDL's variance()
        #
        var = ydata.var()*(float(nx)/float(nx-1))
        if var == 0:
            var = 1.0
        invvar = np.ones(ydata.shape,dtype=ydata.dtype)/var
    if x2 is not None:
        if x2.size != nx:
            raise ValueError('Dimensions of xdata and x2 do not agree.')
    yfit = np.zeros(ydata.shape)
    if invvar.size == 1:
        outmask = True
    else:
        outmask = np.ones(invvar.shape,dtype='bool')
    xsort = xdata.argsort()
    maskwork = (outmask & (invvar > 0))[xsort]
    if 'oldset' in kwargs:
        sset = kwargs['oldset']
        sset.mask = True
        sset.coeff = 0
    else:
        if not maskwork.any():
            raise ValueError('No valid data points.')
            #return (None,None)
        if 'fullbkpt' in kwargs:
            fullbkpt = kwargs['fullbkpt']
        else:
            sset = bspline(xdata[xsort[maskwork]],**kwargs)
            if maskwork.sum() < sset.nord:
                print('Number of good data points fewer than nord.')
                return (sset,outmask)
            if x2 is not None:
                if 'xmin' in kwargs:
                    xmin = kwargs['xmin']
                else:
                    xmin = x2.min()
                if 'xmax' in kwargs:
                    xmax = kwargs['xmax']
                else:
                    xmax = x2.max()
                if xmin == xmax:
                    xmax = xmin + 1
                sset.xmin = xmin
                sset.xmax = xmax
                if 'funcname' in kwargs:
                    sset.funcname = kwargs['funcname']
    xwork = xdata[xsort]
    ywork = ydata[xsort]
    invwork = invvar[xsort]
    if x2 is not None:
        x2work = x2[xsort]
    else:
        x2work = None
    iiter = 0
    error = 0
    qdone = -1
    while (error != 0 or qdone == -1) and iiter <= maxiter:
        goodbk = sset.mask.nonzero()[0]
        if maskwork.sum() <= 1 or not sset.mask.any():
            sset.coeff = 0
            iiter = maxiter + 1
        else:
            if 'requiren' in kwargs:
                i = 0
                while xwork[i] < sset.breakpoints[goodbk[sset.nord]] and i < nx-1:
                    i +=1
                ct = 0
                for ileft in range(sset.nord,sset.mask.sum()-sset.nord+1):
                    while (xwork[i] >= sset.breakpoints[goodbk[ileft]] and
                        xwork[i] < sset.breakpoints[goodbk[ileft+1]] and
                        i < nx-1):
                        ct += invwork[i]*maskwork[i] > 0
                        i += 1
                    if ct >= kwargs['requiren']:
                        ct = 0
                    else:
                        sset.mask[goodbk[ileft]] = False
            error,yfit = sset.fit(xwork,ywork, invwork*maskwork,
                x2=x2work)
        iiter += 1
        inmask = maskwork
        if error == -2:
            return (sset,outmask)
        elif error == 0:
            maskwork,qdone = djs_reject(ywork, yfit, invvar=invwork, inmask=inmask,
                outmask=maskwork, upper=upper, lower=lower)
        else:
            pass
    outmask[xsort] = maskwork
    temp = yfit
    yfit[xsort] = temp
    return (sset,outmask)
