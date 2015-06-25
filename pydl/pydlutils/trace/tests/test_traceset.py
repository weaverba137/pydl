# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from .. import TraceSet, traceset2xy, xy2traceset
from os.path import dirname,join
from astropy.io import fits
#
class TestTraceSet(object):
    def setup(self):
        self.sdss = fits.open(join(dirname(__file__),'t','sdss_traceset.fits'))
        self.boss = fits.open(join(dirname(__file__),'t','boss_traceset.fits'))
        return
    def teardown(self):
        self.sdss.close()
        self.boss.close()
        return
    def test_traceset_sdss(self):
        tset = TraceSet(self.sdss[1])
        assert tset.func == 'legendre'
        assert tset.coeff.shape == (320,5)
    def test_traceset_boss(self):
        tset = TraceSet(self.boss[1])
        assert tset.func == 'legendre'
        assert tset.coeff.shape == (500,6)
