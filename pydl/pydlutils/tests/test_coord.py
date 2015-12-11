# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS
# from astropy.tests.helper import raises
from ..coord import SDSSMuNu, current_mjd, stripe_to_eta, stripe_to_incl


class TestCoord(object):
    """Test the functions in pydl.pydlutils.coord.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_current_mjd(self):
        assert current_mjd() > 50000.0

    def test_munu_to_radec(self):
        munu = SDSSMuNu(mu=0.0*u.deg, nu=0.0*u.deg, stripe=10)
        radec = munu.transform_to(ICRS)
        assert radec.ra.value == 0.0
        assert radec.dec.value == 0.0

    def test_radec_to_munu(self):
        radec = ICRS(ra=0.0*u.deg, dec=0.0*u.deg)
        munu = radec.transform_to(SDSSMuNu(stripe=10))
        assert munu.mu.value == 0.0
        assert munu.nu.value == 0.0
        assert munu.incl.value == 0.0

    def test_stripe_to_eta(self):
        eta = stripe_to_eta(82)
        assert eta == -32.5

    def test_stripe_to_incl(self):
        incl = stripe_to_incl(82)
        assert incl == 0.0
        incl = stripe_to_incl(10)
        assert incl == 0.0
