# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits
from astropy.tests.helper import raises
from astropy.utils.data import get_pkg_data_filename
from ..trace import (fchebyshev, fchebyshev_split, fpoly, func_fit,
                    TraceSet, traceset2xy, xy2traceset)
from .. import PydlutilsException


class TestTrace(object):
    """Test the functions in pydl.pydlutils.trace.
    """

    def setup(self):
        # extracted from spFrame-b1-00057618.fits

        self.sdss = fits.open(get_pkg_data_filename('t/sdss_traceset.fits'))
        # extracted from spFrame-r1-00180406.fits
        self.boss = fits.open(get_pkg_data_filename('t/boss_traceset.fits'))
        return

    def teardown(self):
        self.sdss.close()
        self.boss.close()
        return

    def test_fchebyshev(self):
        x = np.array([-1, -0.5, 0, 0.5, 1], dtype='d')
        #
        # Test order
        #
        with raises(ValueError):
            f = fchebyshev(x, 0)
        #
        # m = 1
        #
        f = fchebyshev(x, 1)
        assert (f == np.ones((1, x.size), dtype='d')).all()
        #
        # m = 2
        #
        f = fchebyshev(x, 2)
        foo = np.ones((2, x.size), dtype='d')
        foo[1, :] = x
        assert np.allclose(f, foo)
        #
        # m = 3
        #
        f = fchebyshev(x, 3)
        foo = np.ones((3, x.size), dtype='d')
        foo[1, :] = x
        foo[2, :] = (2.0*x**2 - 1.0)
        assert np.allclose(f, foo)
        #
        # m = 4
        #
        f = fchebyshev(x, 4)
        foo = np.ones((4, x.size), dtype='d')
        foo[1, :] = x
        foo[2, :] = (2.0*x**2 - 1.0)
        foo[3, :] = (4.0*x**3 - 3.0*x)
        assert np.allclose(f, foo)
        #
        # random float
        #
        f = fchebyshev(2.88, 3)
        assert np.allclose(f, np.array([[1.00], [2.88], [15.5888]]))

    def test_fchebyshev_split(self):
        x = np.array([-1, -0.5, 0, 0.5, 1], dtype='d')
        #
        # Test order
        #
        with raises(ValueError):
            f = fchebyshev_split(x, 0)
        with raises(ValueError):
            f = fchebyshev_split(x, 1)
        #
        # m = 2
        #
        f = fchebyshev_split(x, 2)
        foo = np.ones((2, x.size), dtype='d')
        foo[0, :] = (x >= 0).astype(x.dtype)
        assert np.allclose(f, foo)
        #
        # m = 3
        #
        f = fchebyshev_split(x, 3)
        foo = np.ones((3, x.size), dtype='d')
        foo[0, :] = (x >= 0).astype(x.dtype)
        foo[2, :] = x
        assert np.allclose(f, foo)
        #
        # m = 4
        #
        f = fchebyshev_split(x, 4)
        foo = np.ones((4, x.size), dtype='d')
        foo[0, :] = (x >= 0).astype(x.dtype)
        foo[2, :] = x
        foo[3, :] = (2.0*x**2 - 1.0)
        assert np.allclose(f, foo)
        #
        # m = 5
        #
        f = fchebyshev_split(x, 5)
        foo = np.ones((5, x.size), dtype='d')
        foo[0, :] = (x >= 0).astype(x.dtype)
        foo[2, :] = x
        foo[3, :] = (2.0*x**2 - 1.0)
        foo[4, :] = (4.0*x**3 - 3.0*x)
        assert np.allclose(f, foo)
        #
        # random float
        #
        f = fchebyshev_split(2.88, 3)
        assert np.allclose(f, np.array([[1.00], [1.00], [2.88]]))

    def test_fpoly(self):
        x = np.array([-1, -0.5, 0, 0.5, 1], dtype='d')
        #
        # Test order
        #
        with raises(ValueError):
            f = fpoly(x, 0)
        #
        # m = 1
        #
        f = fpoly(x, 1)
        assert (f == np.ones((1, x.size), dtype='d')).all()
        #
        # m = 2
        #
        f = fpoly(x, 2)
        foo = np.ones((2, x.size), dtype='d')
        foo[1, :] = x
        assert np.allclose(f, foo)
        #
        # m = 3
        #
        f = fpoly(x, 3)
        foo = np.ones((3, x.size), dtype='d')
        foo[1, :] = x
        foo[2, :] = x**2
        assert np.allclose(f, foo)
        #
        # m = 4
        #
        f = fpoly(x, 4)
        foo = np.ones((4, x.size), dtype='d')
        foo[1, :] = x
        foo[2, :] = x**2
        foo[3, :] = x**3
        assert np.allclose(f, foo)
        #
        # random float
        #
        f = fpoly(2.88, 3)
        assert np.allclose(f, np.array([[1.00], [2.88], [8.29440]]))

    def test_func_fit(self):
        np.random.seed(137)
        x = np.linspace(-5, 5, 50)
        y = x**2 + 2*x + 1 + 0.05*np.random.randn(50)
        #
        # Bad inputs
        #
        with raises(ValueError):
            foo = func_fit(x, np.array([1, 2, 3]), 3)
        with raises(ValueError):
            foo = func_fit(x, y, 3, invvar=np.array([0, 0, 0]))
        with raises(ValueError):
            foo = func_fit(x, y, 3, inputfunc=np.array([1.0, 1.0, 1.0]))
        with raises(KeyError):
            foo = func_fit(x, y, 3, function_name='npoly')
        #
        # No good points
        #
        invvar = np.zeros(x.shape, dtype=x.dtype)
        res, yfit = func_fit(x, y, 3, invvar=invvar)
        assert (res == np.zeros((3,), dtype=x.dtype)).all()
        assert (yfit == 0*x).all()
        #
        # One good point
        #
        invvar[2] = 1.0
        res, yfit = func_fit(x, y, 3, invvar=invvar)
        assert (invvar > 0).nonzero()[0][0] == 2
        assert res[0] == y[2]
        assert (yfit == y[2]).all()
        #
        # Various points
        #
        invvar = 1.0/(np.random.random(x.shape)**2)
        assert (invvar < 2).any()
        invvar[invvar < 2] = 0
        res, yfit = func_fit(x, y, 3, invvar=invvar, function_name='poly')
        # assert np.allclose(res,np.array([0.99665423, 1.9945388, 1.00172303]))
        assert np.allclose(res, np.array([0.99996197, 1.99340315, 1.00148004]))
        #
        # Fixed points
        #
        res, yfit = func_fit(x, y, 3, invvar=invvar, function_name='poly',
                            ia=np.array([False, True, True]),
                            inputans=np.array([1.0, 0, 0]))
        # assert np.allclose(res,np.array([1., 1.99454782, 1.00149949]))
        assert np.allclose(res, np.array([1., 1.99340359, 1.00147743]))
        res, yfit = func_fit(x, y, 3, invvar=invvar, function_name='poly',
                            ia=np.array([False, True, False]),
                            inputans=np.array([1.0, 0, 1.0]))
        # assert np.allclose(res,np.array([1., 1.99403239, 1.]))
        assert np.allclose(res, np.array([1., 1.99735654, 1.]))
        #
        # inputfunc
        #
        res, yfit = func_fit(x, y, 3, invvar=invvar, function_name='poly',
                            inputfunc=np.ones(x.shape, dtype=x.dtype))
        # assert np.allclose(res,np.array([0.99665423, 1.9945388, 1.00172303]))
        assert np.allclose(res, np.array([0.99996197, 1.99340315, 1.00148004]))
        #
        # Generate inputans
        #
        y = x**2 + 2*x + 0.05*np.random.randn(50)
        res, yfit = func_fit(x, y, 3, invvar=invvar, function_name='poly',
                            ia=np.array([False, True, True]))
        assert np.allclose(res, np.array([0., 1.99994188, 0.99915111]))

    def test_traceset_sdss(self):
        tset = TraceSet(self.sdss[1].data)
        assert tset.func == 'legendre'
        assert tset.coeff.shape == (320, 5)
        assert not tset.has_jump
        assert tset.xRange == tset.xmax - tset.xmin
        assert tset.nx == int(tset.xmax - tset.xmin + 1)
        assert tset.xmid == 0.5 * (tset.xmin + tset.xmax)
        x, y = traceset2xy(tset)
        tset2 = xy2traceset(x, y, ncoeff=tset.ncoeff)
        assert tset2.xmin == tset.xmin
        assert tset2.xmax == tset.xmax
        assert tset2.func == tset.func
        assert np.allclose(tset2.coeff, tset.coeff)

    def test_traceset_boss(self):
        tset = TraceSet(self.boss[1].data)
        assert tset.func == 'legendre'
        assert tset.coeff.shape == (500, 6)
        assert tset.has_jump
        assert tset.xRange == tset.xmax - tset.xmin
        assert tset.nx == int(tset.xmax - tset.xmin + 1)
        assert tset.xmid == 0.5 * (tset.xmin + tset.xmax)
        x, y = traceset2xy(tset)
        tset2 = xy2traceset(x, y, ncoeff=tset.ncoeff,
                            xjumplo=tset.xjumplo,
                            xjumphi=tset.xjumphi,
                            xjumpval=tset.xjumpval)
        assert tset2.xmin == tset.xmin
        assert tset2.xmax == tset.xmax
        assert tset2.func == tset.func
        assert np.allclose(tset2.coeff, tset.coeff)

    def test_traceset_keywords(self):
        tset = TraceSet(self.boss[1].data)
        x, y = traceset2xy(tset)
        iv = 2.0 * np.ones(x.shape, dtype=x.dtype)
        im = np.ones(x.shape, dtype=np.bool)
        im[1] = False
        tset2 = TraceSet(x, y, func='chebyshev', invvar=iv, maxiter=20,
                         inmask=im)
        assert tset2.func == 'chebyshev'

    def test_traceset_bad(self):
        with raises(PydlutilsException):
            tset = TraceSet(1, 2, 3)
