#
#
#
def skymask(invvar,andmask,ormask=None,ngrow=2):
    """Mask regions where sky-subtraction errors are expected to dominate.
    """
    import numpy as np
    import pydlutils.sdss
    from pydl import smooth
    nrows,npix = invvar.shape
    badmask = np.zeros(invvar.shape,dtype='i4')
    badskychi = pydlutils.sdss.sdss_flagval('SPPIXMASK','BADSKYCHI')
    redmonster = pydlutils.sdss.sdss_flagval('SPPIXMASK','REDMONSTER')
    # brightsky = pydlutils.sdss.sdss_flagval('SPPIXMASK','BRIGHTSKY')
    if ormask is not None:
        badmask = badmask | ((ormask & badskychi) != 0)
        badmask = badmask | ((ormask & redmonster) != 0)
        # badmask = badmask | ((andmask & brightsky) != 0)
    if ngrow > 0:
        width = 2*ngrow + 1
        for k in range(nrows):
            badmask[k,:] = smooth(badmask[k,:]*width, width, True) > 0
    return invvar * (1 - badmask)
