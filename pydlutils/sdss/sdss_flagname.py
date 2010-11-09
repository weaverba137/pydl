#
#
#
def sdss_flagname(flagprefix, flagvalue, concat=False):
    """Return a list of flag names corresponding to the values.

    If concat is True, a single string joined by whitespace will be returned.
    """
    import pydlutils.sdss
    pydlutils.sdss.set_maskbits()
    bitpos = 0
    retval = list()
    while flagvalue > 0:
        bit = flagvalue % 2
        if bit > 0:
            f = filter(lambda x: x[1] == bitpos,pydlutils.sdss.maskbits[flagprefix].items())
            if len(f) > 0:
                retval.append(f[0][0])
        flagvalue >>= 1
        bitpos += 1
    if concat:
        retval = ' '.join(retval)
    return retval

