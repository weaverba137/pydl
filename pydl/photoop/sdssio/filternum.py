# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def filternum(filt='foo'):
    """Return index number for SDSS filters either from a number or name.

    Parameters
    ----------
    filt : str
        The filter name.

    Returns
    -------
    filternum : int
        The corresponding filter number

    Examples
    --------
    >>> filternum('g')
    1
    """
    if filt == 'foo':
        return list(range(5))
    else:
        filters = { 'u':0, 'g':1, 'r':2, 'i':3, 'z':4 }
        return filters[filt]
