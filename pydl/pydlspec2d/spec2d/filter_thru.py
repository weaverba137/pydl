# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def filter_thru(flux,waveimg=None,wset=None,mask=None,filter_prefix='sdss_jun2001',toair=False):
    """Compute throughput in SDSS filters.

    Parameters
    ----------
    flux : array-like
        Spectral flux.
    waveimg : array-like, optional
        Full wavelength solution, with the same shape as `flux`.
    wset : :class:`TraceSet`, optional
        A trace set containing the wavelength solution.  Must be specified
        if `waveimg` is not specified.
    mask : array-like, optional
        Interpolate over pixels where `mask` is non-zero.
    filter_prefix : str, optional
        Specifies a set of filter curves.
    toair : bool, optional
        If ``True``, convert the wavelengths to air from vacuum before computing.

    Returns
    -------
    filter_thru : array-like
        Integrated flux in the filter bands.
    """
    from os import getenv
    from os.path import dirname,join
    import numpy as np
    from astropy.io import ascii
    from astropy.utils.data import get_pkg_data_filename
    from ...goddard.astro import vactoair
    from ...pydlutils.image import djs_maskinterp
    from ...pydlutils.trace import traceset2xy, xy2traceset
    nTrace,nx = flux.shape
    if filter_prefix != 'sdss_jun2001':
        raise ValueError("Filters other than {0} are not available!".format('sdss_jun2001'))
    ffiles = [(join(dirname(__file__),'..','..','pydlutils','data','filters','{0}_{1}_atm.dat'.format(filter_prefix,f))) for f in 'ugriz']
    if waveimg is None and wset is None:
        raise ValueError("Either waveimg or wset must be specified!")
    if waveimg is None:
        pixnorm, logwave = traceset2xy(wset)
        waveimg = 10**logwave
    if toair:
        newwaveimg = vactoair(waveimg)
    else:
        newwaveimg = waveimg
    logwave = np.log10(newwaveimg)
    diffx = np.outer(np.ones((nTrace,),dtype=flux.dtype),np.arange(nx-1,dtype=flux.dtype))
    diffy = logwave[:,1:] - logwave[:,0:nx-1]
    diffset = xy2traceset(diffx,diffy,ncoeff=4,xmin=0,xmax=nx-1)
    pixnorm, logdiff = traceset2xy(diffset)
    logdiff = np.absolute(logdiff)
    if mask is not None:
        flux_interp = djs_maskinterp(flux, mask, iaxis=0)
    res = np.zeros((nTrace,len(ffiles)),dtype=flux.dtype)
    for i,f in enumerate(ffiles):
        filter_data = ascii.read(f,comment='#.*',
            names=('lam', 'respt', 'resbig', 'resnoa', 'xatm'))
        filtimg = logdiff * np.interp(newwaveimg.flatten(),filter_data['lam'].data,filter_data['respt'].data).reshape(logdiff.shape)
        if mask is not None:
            res[:,i] = (flux_interp * filtimg).sum(1)
        else:
            res[:,i] = (flux * filtimg).sum(1)
        sumfilt = filtimg.sum(1)
        res[:,i] = res[:,i] / (sumfilt + (sumfilt <= 0).astype(sumfilt.dtype))
    return res
