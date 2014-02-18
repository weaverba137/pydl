# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_spherematch():
    import numpy as np
    from .. import spherematch
    i1_should_be = np.array([17,  0,  2, 16, 12, 13,  1,  5, 15,  7, 19,  8, 11, 10, 14, 18,  3,
                              9,  6,  4])
    i2_should_be = np.array([ 2,  0, 17,  3, 16, 15, 14,  5,  6, 10,  8, 19, 18,  4,  9,  7, 11,
                              12,  1, 13])
    np.random.seed(137)
    searchrad = 3.0/3600.0
    n = 20
    ra1 = 360.0*np.random.random((n,))
    dec1 = 90.0 - np.rad2deg(np.arccos(2.0*np.random.random((n,)) - 1.0))
    ra2 = ra1 + np.random.normal(0,1.0/3600.0)
    dec2 = dec1 + np.random.normal(0,1.0/3600.0)
    foo = np.arange(n)
    np.random.shuffle(foo)
    i1, i2, d12 = spherematch(ra1, dec1, ra2[foo], dec2[foo], searchrad, maxmatch=0)
    assert (i1 == i1_should_be).all()
    assert (i2 == i2_should_be).all()

