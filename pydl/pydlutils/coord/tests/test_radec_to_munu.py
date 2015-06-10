# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_radec_to_munu():
    from astropy import units as u
    from astropy.coordinates import ICRS
    from .. import SDSSMuNu
    radec = ICRS(ra=0.0*u.deg,dec=0.0*u.deg)
    munu = radec.transform_to(SDSSMuNu(stripe=10))
    assert munu.mu.value == 0.0
    assert munu.nu.value == 0.0
    assert munu.incl.value == 0.0
