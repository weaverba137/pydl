# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..image import djs_maskinterp1, djs_maskinterp


class TestImage(object):
    """Test the functions in pydl.pydlutils.image.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_djs_maskinterp1(self):
        y = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64)
        # Test all good
        yi = djs_maskinterp1(y, np.zeros(y.shape, dtype=np.int32))
        assert (yi == y).all()
        # Test all bad
        yi = djs_maskinterp1(y, np.ones(y.shape, dtype=np.int32))
        assert (yi == y).all()
        # Test one good value
        yt = np.zeros(y.shape, dtype=y.dtype) + y[0]
        yi = djs_maskinterp1(y, np.arange(len(y)))
        assert (yi == yt).all()
        # Test two bad values
        yi = djs_maskinterp1(y, np.array([0, 1, 0, 1, 0]))
        assert np.allclose(y, yi)
        # Because of the default behavior of np.interp(), the 'const'
        # keyword has no actual effect, because it performs the
        # same operation.
        yt = y.copy()
        yt[4] = 3.0
        yi = djs_maskinterp1(y, np.array([0, 1, 0, 0, 1]))
        assert np.allclose(yt, yi)
        yi = djs_maskinterp1(y, np.array([0, 1, 0, 0, 1]), const=True)
        assert np.allclose(yt, yi)
        yt = y.copy()
        yt[0] = 1.0
        yi = djs_maskinterp1(y, np.array([1, 0, 1, 0, 0]), const=True)
        assert np.allclose(yt, yi)
        #
        # Test other x values
        #
        x = np.array([0.0, 0.5, 1.0, 1.5, 2.0], dtype=y.dtype)
        yi = djs_maskinterp1(y, np.array([0, 1, 0, 1, 0]), xval=x)
        assert np.allclose(y, yi)
        yt = y.copy()
        yt[4] = 3.0
        yi = djs_maskinterp1(y, np.array([0, 1, 0, 0, 1]), xval=x)
        assert np.allclose(yt, yi)
        yi = djs_maskinterp1(y, np.array([0, 1, 0, 0, 1]), xval=x, const=True)
        assert np.allclose(yt, yi)
        yt = y.copy()
        yt[0] = 1.0
        yi = djs_maskinterp1(y, np.array([1, 0, 1, 0, 0]), xval=x, const=True)
        assert np.allclose(yt, yi)

    def test_djs_maskinterp(self):
        y = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float64)
        mask = np.array([0, 1, 0])
        with raises(ValueError):
            yi = djs_maskinterp(y, mask)
        mask = np.array([0, 1, 0, 0, 0])
        x = np.array([0.0, 0.5, 1.0], dtype=y.dtype)
        with raises(ValueError):
            yi = djs_maskinterp(y, mask, xval=x)
        # 1-D case
        yi = djs_maskinterp(y, mask)
        assert np.allclose(y, yi)
        # 2-D case
        x = np.array([0.0, 0.5, 1.0, 1.5, 2.0], dtype=y.dtype)
        x = np.vstack((x, x, x))
        y = np.vstack((y, y, y))
        mask = np.vstack((mask, mask, mask))
        with raises(ValueError):
            yi = djs_maskinterp(y, mask)
        with raises(ValueError):
            yi = djs_maskinterp(y, mask, axis=-1)
        with raises(ValueError):
            yi = djs_maskinterp(y, mask, axis=2)
        yi = djs_maskinterp(y, mask, axis=0)
        assert np.allclose(y, yi)
        yi = djs_maskinterp(y, mask, axis=0, xval=x)
        assert np.allclose(y, yi)
        mask[:, 1] = 0
        mask[1, :] = 1
        yi = djs_maskinterp(y, mask, axis=1)
        assert np.allclose(y, yi)
        yi = djs_maskinterp(y, mask, axis=1, xval=x)
        assert np.allclose(y, yi)
        # 3-D case
        x = np.dstack((x, x, x, x, x, x, x))
        y = np.dstack((y, y, y, y, y, y, y))
        mask = np.dstack((mask, mask, mask, mask, mask, mask, mask))
        mask[:, :, :] = 0
        mask[:, :, 5] = 1
        yi = djs_maskinterp(y, mask, axis=0)
        assert np.allclose(y, yi)
        yi = djs_maskinterp(y, mask, axis=0, xval=x)
        assert np.allclose(y, yi)
        mask[:, :, :] = 0
        mask[:, 3, :] = 1
        yi = djs_maskinterp(y, mask, axis=1)
        assert np.allclose(y, yi)
        yi = djs_maskinterp(y, mask, axis=1, xval=x)
        assert np.allclose(y, yi)
        mask[:, :, :] = 0
        mask[1, :, :] = 1
        yi = djs_maskinterp(y, mask, axis=2)
        assert np.allclose(y, yi)
        yi = djs_maskinterp(y, mask, axis=2, xval=x)
        assert np.allclose(y, yi)
        # 4-D case
        y = np.random.random((2, 2, 2, 2))
        with raises(ValueError):
            yi = djs_maskinterp(y, (y > 0.5), axis=0)
