def aesthetics(flux,invvar,**kwargs):
    import numpy as np
    from scipy.special import erf
    from pydlutils.image import djs_maskinterp
    from pydlspec2d import Pydlspec2dException
    badpts = invvar == 0
    if badpts.any():
        if 'method' not in kwargs:
            method = 'traditional'
        else:
            method = kwargs['method']
        if method == 'traditional':
            newflux = djs_maskinterp(flux,invvar == 0,const=True)
        elif method == 'noconst':
            newflux = djs_maskinterp(flux,invvar == 0)
        elif method == 'mean':
            newflux = flux.copy()
            goodpts = invvar > 0
            newflux[~goodpts] = newflux[goodpts].mean()
        elif method == 'damp':
            l = 250 # damping length in pixels
            goodpts = invvar.nonzero()[0]
            nflux = flux.size
            mingood = goodpts.min()
            maxgood = goodpts.max()
            newflux = djs_maskinterp(flux,invvar == 0,const=True)
            pixels = np.arange(nflux,dtype='f')
            if mingood > 0:
                damp1 = float(min(mingood,l))
                newflux *= 0.5*(1.0+erf((pixels-mingood)/damp1))
            if maxgood < (nflux - 1):
                damp2 = float(min(maxgood,l))
                newflux *= 0.5*(1.0+erf((maxgood-pixels)/damp2))
        elif method == 'nothing':
            newflux = flux.copy()
        else:
            raise Pydlspec2dException("Unknown method: %s" % method)
        return newflux
    else:
        return flux

