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
    import numpy as np
    from ..goddard.astro import vactoair
    from ..pydlutils.trace import traceset2xy, xy2traceset
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
    
    return
