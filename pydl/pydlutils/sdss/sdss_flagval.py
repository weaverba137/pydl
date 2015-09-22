# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.extern import six

def sdss_flagval(flagname,bitname):
    """Convert bitmask names into values.

    Converts human-readable bitmask names into numerical values.  The inputs
    are not case-sensitive; all inputs are converted to upper case internally.

    Parameters
    ----------
    flagname : :class:`str`
        The name of a bitmask group.
    bitname : :class:`str` or :class:`list`
        The name(s) of the specific bitmask(s) within the `flagname` group.

    Returns
    -------
    sdss_flagval : :class:`numpy.uint64`
        The value of the bitmask name(s).

    Raises
    ------
    KeyError
        If `flagname` or `bitname` are invalid names.

    Examples
    --------
    >>> from pydl.pydlutils.sdss import sdss_flagval
    >>> sdss_flagval('ANCILLARY_TARGET1',['BLAZGX','ELG','BRIGHTGAL']) # doctest: +REMOTE_DATA
    2310346608843161600
    """
    from numpy import uint64
    from . import maskbits
    if maskbits is None: # pragma: no cover
        from .set_maskbits import set_maskbits
        maskbits = set_maskbits()
    #
    # Make sure inlabel is a list
    #
    if isinstance(bitname,six.string_types):
        bitnames = [bitname.upper()]
    else:
        bitnames = [b.upper() for b in bitname]
    flagu = flagname.upper()
    flagvalue = uint64(0)
    for bit in bitnames:
        if flagu in maskbits:
            if bit in maskbits[flagu]:
                flagvalue += uint64(2)**uint64(maskbits[flagu][bit])
            else:
                raise KeyError("Unknown bit label {0} for flag group {1}!".format(bit, flagu))
        else:
            raise KeyError("Unknown flag group {0}!".format(flagu))
    return flagvalue
