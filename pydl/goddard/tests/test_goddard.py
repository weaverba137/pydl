# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..astro import airtovac, gcirc, get_juldate, vactoair
from ..math import flegendre
from ..misc import cirrange


class TestGoddard(object):
    """Test the goddard package.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_airtovac(self):
        vacuum = airtovac(1900.0)
        assert vacuum == 1900.0
        vacuum = airtovac(2000.0)
        assert np.allclose(vacuum, 2000.6475)
        air = np.array([1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0])
        vacuum = airtovac(air)
        assert np.allclose(vacuum,
                           np.array([1800.0, 1900.0, 2000.6475,
                                     2100.6664, 2200.6868, 2300.7083]))
        vacuum = airtovac(6056.125)
        assert np.allclose(vacuum, 6057.8019)
        #
        # Regression test for #8.
        #
        wave = air.reshape(2, 3)
        vacuum = airtovac(wave)
        assert np.allclose(vacuum,
                           np.array([[1800.0, 1900.0, 2000.6475],
                                     [2100.6664, 2200.6868, 2300.7083]]))

    def test_cirrange(self):
        ra1 = np.linspace(-4.0*np.pi, 4.0*np.pi, 100)
        ra2 = cirrange(ra1, radians=True)
        assert (ra2 == (ra1 % (2.0 * np.pi))).all()
        ra1 = np.rad2deg(ra1)
        ra2 = cirrange(ra1)
        assert (ra2 == (ra1 % 360.0)).all()

    def test_flegendre(self):
        x = np.array([-1, -0.5, 0, 0.5, 1], dtype='d')
        #
        # Test order
        #
        with raises(ValueError):
            f = flegendre(x, 0)
        #
        # m = 1
        #
        f = flegendre(x, 1)
        assert (f == np.ones((1, x.size), dtype='d')).all()
        #
        # m = 2
        #
        f = flegendre(x, 2)
        foo = np.ones((2, x.size), dtype='d')
        foo[1, :] = x
        assert np.allclose(f, foo)
        #
        # m = 3
        #
        f = flegendre(x, 3)
        foo = np.ones((3, x.size), dtype='d')
        foo[1, :] = x
        foo[2, :] = 0.5 * (3.0*x**2 - 1.0)
        assert np.allclose(f, foo)
        #
        # m = 4
        #
        f = flegendre(x, 4)
        foo = np.ones((4, x.size), dtype='d')
        foo[1, :] = x
        foo[2, :] = 0.5*(3.0*x**2 - 1.0)
        foo[3, :] = 0.5*(5.0*x**3 - 3.0*x)
        assert np.allclose(f, foo)
        #
        # random float
        #
        f = flegendre(2.88, 3)
        assert np.allclose(f, np.array([[1.00], [2.88], [11.9416]]))

    def test_gcirc(self):
        np.random.seed(137)
        #
        # Start in radians
        #
        offset = 5.0e-6  # approx 1 arcsec
        ra1 = 2.0 * np.pi * np.random.rand(100)
        dec1 = np.pi/2.0 - np.arccos(2.0*np.random.rand(100) - 1.0)
        ra2 = ra1 + offset
        ra2 = np.where((ra2 > 2.0*np.pi), ra2 - 2.0*np.pi, ra2)
        dec2 = np.where((dec1 > 0), dec1 - offset, dec1 + offset)
        deldec2 = (dec2-dec1)/2.0
        delra2 = (ra2-ra1)/2.0
        sindis = np.sqrt(np.sin(deldec2) * np.sin(deldec2) +
                         np.cos(dec1) * np.cos(dec2) *
                         np.sin(delra2) * np.sin(delra2))
        dis = 2.0*np.arcsin(sindis)
        #
        # units = 0
        #
        d0 = gcirc(ra1, dec1, ra2, dec2, units=0)
        assert np.allclose(d0, dis)
        #
        # units = 2
        #
        d0 = gcirc(np.rad2deg(ra1)/15.0, np.rad2deg(dec1),
                   np.rad2deg(ra2)/15.0, np.rad2deg(dec2), units=1)
        assert np.allclose(d0, np.rad2deg(dis)*3600.0)
        #
        # units = 2
        #
        d0 = gcirc(np.rad2deg(ra1), np.rad2deg(dec1), np.rad2deg(ra2),
                   np.rad2deg(dec2), units=2)
        assert np.allclose(d0, np.rad2deg(dis)*3600.0)
        #
        # Units = whatever
        #
        with raises(ValueError):
            d0 = gcirc(ra1, dec1, ra2, dec2, units=5)

    def test_get_juldate(self):
        now = get_juldate()
        assert now > 2400000.5
        assert get_juldate(0) == 40587.0 + 2400000.5
        assert get_juldate(86400) == 40588.0 + 2400000.5

    def test_vactoair(self):
        air = vactoair(1900.0)
        assert air == 1900.0
        air = vactoair(2000.0)
        assert np.allclose(air, 1999.3526)
        vacuum = np.array([1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0])
        air = vactoair(vacuum)
        assert np.allclose(air, np.array([1800.0, 1900.0, 1999.3526,
                                          2099.3337, 2199.3133, 2299.2918]))
        #
        # Regression test for #8.
        #
        wave = vacuum.reshape(2, 3)
        air = vactoair(wave)
        assert np.allclose(air,
                           np.array([[1800.0, 1900.0, 1999.3526],
                                     [2099.3337, 2199.3133, 2299.2918]]))
