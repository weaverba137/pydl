# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_flagexist(flagname, bitname, flagexist=False):
    """Check for the existence of flags.

    Parameters
    ----------
    flagname : str
        The name of a bitmask group. Not case-sensitive.
    bitname : str or list
        The name(s) of the specific bitmask(s) within the `flagname` group.
    flagexist : bool
        If flagexist is True, return a tuple with the second component indicating
        whether the binary flag named `flagname` exists, even if `bitname` is wrong.

    Returns
    -------
    sdss_flagexist : bool or tuple
        A boolean value or a tuple of bool.
    """
    from . import maskbits
    #
    # Make sure label is a list
    #
    #if isinstance(bitname,str):
    #    bitnames = [bitname.upper()]
    #else:
    #    bitnames = [b.upper() for b in bitname]
    f = False
    l = False
    if flagname.upper() in maskbits:
        f = True
        if bitname in maskbits[flagname.upper()]:
            l = True
    if flagexist:
        return (l,f)
    else:
        return l
