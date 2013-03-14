# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def skymask(invvar,andmask,ormask=None,ngrow=2):
    """Mask regions where sky-subtraction errors are expected to dominate.
    """
    from numpy import zeros
    from ...pydlutils.sdss import sdss_flagval
    from ... import smooth
    nrows,npix = invvar.shape
    badmask = zeros(invvar.shape,dtype='i4')
    badskychi = sdss_flagval('SPPIXMASK','BADSKYCHI')
    redmonster = sdss_flagval('SPPIXMASK','REDMONSTER')
    # brightsky = sdss_flagval('SPPIXMASK','BRIGHTSKY')
    if ormask is not None:
        badmask = badmask | ((ormask & badskychi) != 0)
        badmask = badmask | ((ormask & redmonster) != 0)
        # badmask = badmask | ((andmask & brightsky) != 0)
    if ngrow > 0:
        width = 2*ngrow + 1
        for k in range(nrows):
            badmask[k,:] = smooth(badmask[k,:]*width, width, True) > 0
    return invvar * (1 - badmask)
