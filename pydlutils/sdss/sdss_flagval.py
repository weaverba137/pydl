#
#
#
def sdss_flagval(flagprefix,inlabel):
    """Return bitmask values corresponding to labels.

    Some bitfields are larger than 32 bits, so use python long() just to be safe.
    """
    pydlutils.sdss.set_maskbits()
    #
    # Make sure inlabel is a list
    #
    if isinstance(inlabel,str):
        inlabel = [inlabel]
    flagvalue = 0L
    for l in inlabel:
        if flagprefix.upper() in maskbits:
            if l.upper() in maskbits[flagprefix.upper()]:
                flagvalue += 2L**maskbits[flagprefix.upper()][l.upper()]
            else:
                print "Unknown bit label %s for flag %s!" % (flagprefix.upper(), l.upper())
                return None
        else:
            print "Unknown flag %s!" % flagprefix.upper()
            return None
    return flagvalue
