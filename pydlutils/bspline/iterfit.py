#
#
#
def iterfit(xdata,ydata,invvar=None,**kwargs):
    """Iteratively fit a bspline set to data.
    """
    import numpy as np
    from pydlutils.bspline import bspline
    from pydlutils.math import djs_reject
    nx = xdata.size
    if ydata.size != nx:
        raise ValueError('Dimensions of xdata and ydata do not agree.')
    if 'upper' in kwargs:
        upper = kwargs['upper']
    else:
        upper = 5
    if 'lower' in kwargs:
        lower = kwargs['lower']
    else:
        lower = 5
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
    if 'x2' in kwargs:
        if kwargs['x2'].size != nx:
            raise ValueError('Dimensions of xdata and x2 do not agree.')
    if 'maxiter' in kwargs:
        maxiter = kwargs['maxiter']
    else:
        maxiter = 10
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
            print 'No valid data points.'
            return (None,None)
        if 'fullbkpt' in kwargs:
            fullbkpt = kwargs['fullbkpt']
        else:
            sset = bspline(xdata[xsort[maskwork]],**kwargs)
            if maskwork.sum() < sset.nord:
                print 'Number of good data points fewer than nord.'
                return (sset,outmask)
            if 'x2' in kwargs:
                if 'xmin' in kwargs:
                    xmin = kwargs['xmin']
                else:
                    xmin = kwargs['x2'].min()
                if 'xmax' in kwargs:
                    xmax = kwargs['xmax']
                else:
                    xmax = kwargs['x2'].max()
                if xmin == xmax:
                    xmax = xmin + 1
                sset.xmin = xmin
                sset.xmax = xmax
                if 'funcname' in kwargs:
                    sset.funcname = kwargs['funcname']
    xwork = xdata[xsort]
    ywork = ydata[xsort]
    invwork = invvar[xsort]
    if 'x2' in kwargs:
        x2work = kwargs['x2'][xsort]
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
                x2=x2work, **kwargs)
        iiter += 1
        inmask = maskwork
        if error == -2:
            return (sset,outmask)
        elif error == 0:
            maskwork,qdone = djs_reject(ywork, yfit, invvar=invwork, inmask=inmask,
                outmask=maskwork, upper=upper, lower=lower, **kwargs)
        else:
            pass
    outmask[xsort] = maskwork
    temp = yfit
    yfit[xsort] = temp
    return (sset,outmask)

