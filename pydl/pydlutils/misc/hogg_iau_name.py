# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def hogg_iau_name(ra,dec,prefix='SDSS',precision=1):
    """Properly format astronomical source names to the IAU convention.

    Parameters
    ----------
    ra : :class:`float` or :class:`numpy.ndarray`
        Right ascencion in decimal degrees
    dec : :class:`float` or :class:`numpy.ndarray`
        Declination in decimal degrees.
    prefix : :class:`str`, optional
        Add this prefix to the string, defaults to 'SDSS'.
    precision : :class:`int`, optional
        Display this many digits of precision on seconds, default 1.

    Returns
    -------
    hogg_iau_name : :class:`str` or :class:`list`
        The IAU name for the coordinates.

    Examples
    --------
    >>> from pydl.pydlutils.misc import hogg_iau_name
    >>> hogg_iau_name(354.120375,-0.544777778)
    'SDSS J233628.89-003241.2'
    """
    import numpy as np
    #
    # Promote scalar values to arrays.
    #
    if isinstance(ra,float):
        ra = np.array([ra])
    if isinstance(dec,float):
        dec = np.array([dec])
    h = ra/15.0
    rah = np.floor(h)
    ram = np.floor(60.0*(h-rah))
    ras = 60.0*(60.0*(h-rah) - ram)
    ras = np.floor(ras*10.0**(precision+1))/10.0**(precision+1)
    rasformat = "{{2:0{0:d}.{1:d}f}}".format(precision+4, precision+1)
    rah = rah.astype(np.int32)
    ram = ram.astype(np.int32)
    desgn = np.array(list('+'*len(dec)))
    desgn[dec < 0] = '-'
    adec = np.absolute(dec)
    ded = np.floor(adec)
    dem = np.floor(60.0*(adec-ded))
    des = 60.0*(60.0*(adec-ded) - dem)
    des = np.floor(des*10.0**precision)/10.0**precision
    desformat = "{{6:0{0:d}.{1:d}f}}".format(precision+3, precision)
    if precision == 0:
        desformat = "{6:02d}"
        des = des.astype(np.int32)
    ded = ded.astype(np.int32)
    dem = dem.astype(np.int32)
    adformat = "{{0:02d}}{{1:02d}}{ras}{{3:s}}{{4:02d}}{{5:02d}}{des}".format(ras=rasformat,des=desformat)
    #
    # The easy way doesn't work if numpy version is less than 1.7.0.  Prior to this
    # version, numpy scalars didn't support __format__.  See
    # http://projects.scipy.org/numpy/ticket/1675
    #
    try:
        adstr = [adformat.format(*x) for x in zip(rah, ram, ras, desgn, ded, dem, des)]
    except ValueError:
        adstr = [adformat.format(*x) for x in zip(
            rah.tolist(), ram.tolist(), ras.tolist(), desgn.tolist(),
            ded.tolist(), dem.tolist(), des.tolist())]
    if prefix == '':
        jstr = 'J'
    else:
        jstr = ' J'
    name = ["{0}{1}{2}".format(prefix, jstr, x) for x in adstr]
    if len(ra) == 1:
        return name[0]
    else:
        return name
#
#
#
def hogg_iau_name_main(): # pragma: no cover
    from astropy.utils.compat import argparse
    parser = argparse.ArgumentParser(description='Properly format astronomical source names to the IAU convention.')
    parser.add_argument('-P', '--precision', dest='precision', action='store',
        metavar='N', default=1, type=int, help='Digits of precision to add to the declination.')
    parser.add_argument('-p', '--prefix', dest='prefix', action='store',
        metavar='STR', default='SDSS', help='Add this prefix to the name.')
    parser.add_argument('ra', metavar='RA', type=float,
        help='Right Ascension.')
    parser.add_argument('dec', metavar='Dec', type=float,
        help='Declination.')
    options = parser.parse_args()
    print(hogg_iau_name(options.ra,options.dec,
        prefix=options.prefix,precision=options.precision))
    return 0
