# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_flagexist(flagprefix, label, flagexist=False):
    """Returns true if a flag exists. If flagexist is True,
    return a tuple with the second component indicating
    whether the binary flag exists, even if the label is wrong.
    """
    from . import maskbits
    #set_maskbits()
    #
    # Make sure label is a list
    #
    #if isinstance(label,str):
    #    label = [label]
    f = False
    l = False
    if flagprefix in maskbits:
        f = True
        if label in maskbits[flagprefix]:
            l = True
    if flagexist:
        return (l,f)
    else:
        return l
