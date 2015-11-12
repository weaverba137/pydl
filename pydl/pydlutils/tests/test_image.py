# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from ..image import djs_maskinterp1


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
