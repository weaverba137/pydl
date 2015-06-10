# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def wavevector(minfullwave,maxfullwave,zeropoint=3.5,binsz=1.0e-4,wavemin=None,wavemax=None):
    """Return an array of wavelengths.

    Parameters
    ----------
    minfullwave : float
        Minimum wavelength.
    maxfullwave : float
        Maximum wavelength.
    zeropoint : float, optional
        Offset of the input wavelength values.
    binsz : float, optional
        Separation between wavelength values.
    wavemin : float, optional
        If this is set the values of `minfullwave` and `zeropoint` are ignored.
    wavemax : float, optional
        If this is set the value of `maxfullwave` is ignored.

    Returns
    -------
    wavevector : numpy.ndarray
        Depending on the values of `minfullwave`, `binsz`, etc., the resulting
        array could be interpreted as an array of wavelengths or an array of
        log(wavelength).
    """
    from numpy import arange
    if wavemin is not None:
        spotmin = 0
        if wavemax is not None:
            spotmax = int((wavemax - wavemin)/binsz)
        else:
            spotmax = int((maxfullwave - wavemin)/binsz)
            wavemax = spotmax * binsz + wavemin
    else:
        spotmin = int((minfullwave - zeropoint)/binsz) + 1
        spotmax = int((maxfullwave - zeropoint)/binsz)
        wavemin = spotmin * binsz + zeropoint
        wavemax = spotmax * binsz + zeropoint
    nfinalpix = spotmax - spotmin + 1
    finalwave = arange(nfinalpix,dtype='d') *binsz + wavemin
    return finalwave
