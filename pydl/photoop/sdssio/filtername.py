# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def filtername(f):
    """Return the name of a filter given its number.

    Parameters
    ----------
    f : int
        The filter number.

    Returns
    -------
    filtername : str
        The corresponding filter name.

    Examples
    --------
    >>> filtername(0)
    'u'
    """
    if isinstance(f,str):
        return f
    fname = ('u','g','r','i','z')
    return fname[f]
