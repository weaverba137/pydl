# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_flagexist(flagname, bitname, flagexist=False, whichexist=False):
    """Check for the existence of flags.

    Parameters
    ----------
    flagname : :class:`str`
        The name of a bitmask group. Not case-sensitive.
    bitname : :class:`str` or :class:`list`
        The name(s) of the specific bitmask(s) within the `flagname` group.
    flagexist : :class:`bool`, optional
        If flagexist is True, return a tuple with the second component indicating
        whether the binary flag named `flagname` exists, even if `bitname` is wrong.
    whichexist : :class:`bool`, optional
        If whichexist is True, return a list containing existence test results
        for each individual flag.

    Returns
    -------
    sdss_flagexist : :class:`bool` or :func:`tuple`
        A boolean value or a tuple of bool.
    """
    from . import maskbits
    if maskbits is None: # pragma: no cover
        from .set_maskbits import set_maskbits
        maskbits = set_maskbits()
    #
    # Make sure label is a list
    #
    if isinstance(bitname,str):
        bitnames = [bitname.upper()]
    else:
        bitnames = [b.upper() for b in bitname]
    f = False
    l = False
    which = [False]*len(bitnames)
    if flagname.upper() in maskbits:
        f = True
        which = [n in maskbits[flagname.upper()] for n in bitnames]
        l = sum(which) == len(which)
    if flagexist and whichexist:
        return (l,f,which)
    elif flagexist:
        return (l,f)
    elif whichexist:
        return (l,which)
    else:
        return l
