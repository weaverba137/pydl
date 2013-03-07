#
#
#
def sdss_flagexist(flagprefix, label, flagexist=False):
    """Returns true if a flag exists. If flagexist is True,
    return a tuple with the second component indicating
    whether the binary flag exists, even if the label is wrong.
    """
    import pydlutils.sdss
    pydlutils.sdss.set_maskbits()
    #
    # Make sure label is a list
    #
    #if isinstance(label,str):
    #    label = [label]
    f = False
    l = False
    if flagprefix in pydlutils.sdss.maskbits:
        f = True
        if label in pydlutils.sdss.maskbits[flagprefix]:
            l = True
    if flagexist:
        return (l,f)
    else:
        return l
