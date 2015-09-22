# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def filtername(f):
    """Return the name of a filter given its number.

    Parameters
    ----------
    f : :class:`int`
        The filter number.

    Returns
    -------
    filtername : :class:`str`
        The corresponding filter name.

    Examples
    --------
    >>> filtername(0)
    'u'
    """
    from astropy.extern.six import string_types
    if isinstance(f,string_types):
        return f
    fname = ('u','g','r','i','z')
    return fname[f]
