# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_flagname(flagname, flagvalue, concat=False):
    """Return a list of flag names corresponding to the values.

    Parameters
    ----------
    flagname : str
        The name of a bitmask group. Not case-sensitive.
    flagvalue : long
        The value to be converted into bitmask names.
    concat : bool, optional
        If set to ``True``, the list of names is converted to a space-separated string.

    Returns
    -------
    sdss_flagname : str or list
        The names of the bitmasks encoded in `flagvalue`.

    Raises
    ------
    KeyError
        If `flagname` is invalid

    Examples
    --------
    >>> pydl.pydlutils.sdss.sdss_flagname('ANCILLARY_TARGET1',2310346608843161600L)
    ['BRIGHTGAL', 'BLAZGX', 'ELG']
    """
    from . import maskbits
    flagu = flagname.upper()
    bitpos = 0
    retval = list()
    while flagvalue > 0:
        bit = flagvalue % 2
        if bit > 0:
            try:
                f = filter(lambda x: x[1] == bitpos,maskbits[flagu].items())
            except KeyError:
                raise KeyError("Unknown flag group {0}!".format(flagu))
            if len(f) > 0:
                retval.append(f[0][0])
        flagvalue >>= 1
        bitpos += 1
    if concat:
        retval = ' '.join(retval)
    return retval

