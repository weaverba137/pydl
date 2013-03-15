# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def spec_append(spec1, spec2, pixshift=0):
    """Append the array spec2 to the array spec1 & return a new array.

    If the dimension of these arrays is the same, then append as [spec1,spec2].
    If not, increase the size of the smaller array & fill with zeros.

    Parameters
    ----------
    spec1, spec2 : ndarray
        Append `spec2` to `spec1`.
    pixshift : int, optional
        If `pixshift` is set to a positive integer, `spec2` will be padded with
        `pixshift` zeros on the left size.  If `pixshift` is set to a
        negative integer, `spec1` will be padded with abs(`pixshift`) zeros on
        the left side.  If not set, all zeros will be padded on the right side.

    Returns
    -------
    spec_append : ndarray
        A new array containing both `spec1` and `spec2`.
    """
    from numpy import zeros
    nrows1,npix1 = spec1.shape
    nrows2,npix2 = spec2.shape
    nrows = nrows1+nrows2
    nadd1 = 0
    nadd2 = 0
    if pixshift != 0:
        if pixshift < 0:
            nadd1 = -pixshift
        else:
            nadd2 = pixshift
    maxpix = max(npix1+nadd1,npix2+nadd2)
    spec3 = zeros((nrows,maxpix),dtype=spec1.dtype)
    spec3[0:nrows1,nadd1:nadd1+npix1] = spec1
    spec3[nrows1:nrows,nadd2:nadd2+npix2] = spec2
    return spec3
