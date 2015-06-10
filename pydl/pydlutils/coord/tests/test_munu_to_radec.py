# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_munu_to_radec():
    from astropy import units as u
    from astropy.coordinates import ICRS
    from .. import SDSSMuNu
    munu = SDSSMuNu(mu=0.0*u.deg,nu=0.0*u.deg,stripe=10)
    radec = munu.transform_to(ICRS)
    assert radec.ra.value == 0.0
    assert radec.dec.value == 0.0
