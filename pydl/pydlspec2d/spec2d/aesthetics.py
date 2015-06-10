# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def aesthetics(flux,invvar,method='traditional'):
    """Add nice values to a spectrum where it is masked.

    Parameters
    ----------
    flux : ndarray
        The spectrum to clean up.
    invvar : ndarray
        Inverse variance of the spectrum.
    method : { 'traditional', 'noconst', 'mean', 'damp', 'nothing' }, optional
        Apply this method to clean up the spectrum.  Default is 'traditional'.

    Returns
    -------
    aesthetics : ndarray
        A cleaned-up spectrum.
    """
    import numpy as np
    from scipy.special import erf
    from ...pydlutils.image import djs_maskinterp
    from .. import Pydlspec2dException
    badpts = invvar == 0
    if badpts.any():
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
            raise Pydlspec2dException("Unknown method: {0}".format(method))
        return newflux
    else:
        return flux
