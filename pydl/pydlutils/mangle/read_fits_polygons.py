# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def read_fits_polygons(filename):
    """Read a "polygon" format FITS file.

    Unlike the IDL version of this routine, the ``XCAPS`` and ``CMPCAPS`` columns
    are *not* replaced with ``poly['CAPS']['X']`` and ``poly['CAPS']['CM']``.

    Parameters
    ----------
    filename : :class:`str`
        Name of FITS file to read.

    Returns
    -------
    read_fits_polygons : :class:`~astropy.io.fits.FITS_rec`
        The data contained in HDU 1 of the FITS file.
    """
    from astropy.io import fits
    with fits.open(filename) as hdulist:
        poly = hdulist[1].data
    # poly['CAPS'] = {'CM':list(), 'X':list()}
    # for k in range(len(poly['NCAPS'])):
    #     poly['CAPS']['CM'].append(poly['CMCAPS'][k,0:poly['NCAPS'][k]])
    #     poly['CAPS']['X'].append(poly['XCAPS'][k].reshape((8,3))[0:poly['NCAPS'][k]])
    # del poly['CMCAPS']
    # del poly['XCAPS']
    return poly
