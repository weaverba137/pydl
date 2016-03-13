# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


def wavevector(minfullwave, maxfullwave, zeropoint=3.5, binsz=1.0e-4,
               wavemin=None):
    """Return an array of wavelengths.

    Parameters
    ----------
    minfullwave : :class:`float`
        Minimum wavelength.
    maxfullwave : :class:`float`
        Maximum wavelength.
    zeropoint : :class:`float`, optional
        Offset of the input wavelength values.
    binsz : :class:`float`, optional
        Separation between wavelength values.
    wavemin : :class:`float`, optional
        If this is set the values of `minfullwave` and `zeropoint` are ignored.

    Returns
    -------
    :class:`numpy.ndarray`
        Depending on the values of `minfullwave`, `binsz`, etc., the resulting
        array could be interpreted as an array of wavelengths or an array of
        log(wavelength).
    """
    from numpy import arange
    if wavemin is not None:
        spotmin = 0
        spotmax = int((maxfullwave - wavemin)/binsz)
        wavemax = spotmax * binsz + wavemin
    else:
        spotmin = int((minfullwave - zeropoint)/binsz) + 1
        spotmax = int((maxfullwave - zeropoint)/binsz)
        wavemin = spotmin * binsz + zeropoint
        wavemax = spotmax * binsz + zeropoint
    nfinalpix = spotmax - spotmin + 1
    finalwave = arange(nfinalpix, dtype='d') * binsz + wavemin
    return finalwave
