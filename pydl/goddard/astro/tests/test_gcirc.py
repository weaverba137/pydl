# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_gcirc():
    from .. import gcirc
    import numpy as np
    from astropy.tests.helper import raises
    np.random.seed(137)
    #
    # Start in radians
    #
    offset = 5.0e-6 # approx 1 arcsec
    ra1 = 2.0*np.pi*np.random.rand(100)
    dec1 = np.pi/2.0 - np.arccos(2.0*np.random.rand(100) - 1.0)
    ra2 = ra1 + offset
    ra2 = np.where((ra2 > 2.0*np.pi), ra2 - 2.0*np.pi, ra2)
    dec2 = np.where((dec1 > 0), dec1 - offset, dec1 + offset)
    deldec2 = (dec2-dec1)/2.0
    delra2 =  (ra2-ra1)/2.0
    sindis = np.sqrt( np.sin(deldec2)*np.sin(deldec2) +
        np.cos(dec1)*np.cos(dec2)*np.sin(delra2)*np.sin(delra2) )
    dis = 2.0*np.arcsin(sindis)
    #
    # units = 0
    #
    d0 = gcirc(ra1,dec1,ra2,dec2,units=0)
    assert np.allclose(d0,dis)
    #
    # units = 2
    #
    d0 = gcirc(np.rad2deg(ra1)/15.0,np.rad2deg(dec1),np.rad2deg(ra2)/15.0,np.rad2deg(dec2),units=1)
    assert np.allclose(d0,np.rad2deg(dis)*3600.0)
    #
    # units = 2
    #
    d0 = gcirc(np.rad2deg(ra1),np.rad2deg(dec1),np.rad2deg(ra2),np.rad2deg(dec2),units=2)
    assert np.allclose(d0,np.rad2deg(dis)*3600.0)
    #
    # Units = whatever
    #
    with raises(ValueError):
        d0 = gcirc(ra1,dec1,ra2,dec2,units=5)
