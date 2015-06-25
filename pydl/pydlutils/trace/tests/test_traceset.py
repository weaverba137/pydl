# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from numpy import allclose
from .. import TraceSet, traceset2xy, xy2traceset
from ... import PydlutilsException
from os.path import dirname,join
from astropy.io import fits
from astropy.tests.helper import raises
#
class TestTraceSet(object):
    def setup(self):
        # extracted from spFrame-b1-00057618.fits
        self.sdss = fits.open(join(dirname(__file__),'t','sdss_traceset.fits'))
        # extracted from spFrame-r1-00180406.fits
        self.boss = fits.open(join(dirname(__file__),'t','boss_traceset.fits'))
        return
    def teardown(self):
        self.sdss.close()
        self.boss.close()
        return
    def test_traceset_sdss(self):
        tset = TraceSet(self.sdss[1].data)
        assert tset.func == 'legendre'
        assert tset.coeff.shape == (320,5)
        assert not tset.has_jump
        assert tset.xRange == tset.xmax - tset.xmin
        assert tset.nx == int(tset.xmax - tset.xmin + 1)
        assert tset.xmid == 0.5 * ( tset.xmin + tset.xmax)
        x,y = traceset2xy(tset)
        tset2 = xy2traceset(x,y,ncoeff=tset.ncoeff)
        assert tset2.xmin == tset.xmin
        assert tset2.xmax == tset.xmax
        assert tset2.func == tset.func
        assert allclose(tset2.coeff,tset.coeff)
    def test_traceset_boss(self):
        tset = TraceSet(self.boss[1].data)
        assert tset.func == 'legendre'
        assert tset.coeff.shape == (500,6)
        assert tset.has_jump
        assert tset.xRange == tset.xmax - tset.xmin
        assert tset.nx == int(tset.xmax - tset.xmin + 1)
        assert tset.xmid == 0.5 * ( tset.xmin + tset.xmax)
        x,y = traceset2xy(tset)
        tset2 = xy2traceset(x,y,ncoeff=tset.ncoeff,
            xjumplo=tset.xjumplo,xjumphi=tset.xjumphi,xjumpval=tset.xjumpval)
        assert tset2.xmin == tset.xmin
        assert tset2.xmax == tset.xmax
        assert tset2.func == tset.func
        assert allclose(tset2.coeff,tset.coeff)
    def test_traceset_bad(self):
        with raises(PydlutilsException):
            tset = TraceSet(1,2,3)
