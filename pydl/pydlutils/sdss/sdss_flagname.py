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
    >>> from pydl.pydlutils.sdss import sdss_flagname
    >>> sdss_flagname('ANCILLARY_TARGET1',2310346608843161600)
    ['BRIGHTGAL', 'BLAZGX', 'ELG']
    """
    from . import maskbits
    if maskbits is None:
        from .set_maskbits import set_maskbits
        maskbits = set_maskbits()
    flagu = flagname.upper()
    bits = list()
    for bit in range(64):
        if (flagvalue & (1 << bit)) != 0:
            bits.append(bit)
    retval = list()
    for bit in bits:
        try:
            f = [ x for x in maskbits[flagu].items() if x[1] == bit ]
        except KeyError:
            raise KeyError("Unknown flag group {0}!".format(flagu))
        if f:
            retval.append(f[0][0])
    if concat:
        retval = ' '.join(retval)
    return retval

