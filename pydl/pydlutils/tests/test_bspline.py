# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..bspline import iterfit
from ... import smooth


class TestBspline(object):
    """Test the functions in pydl.pydlutils.bspline.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_iterfit(self):
        y0 = np.array([0.661984, 0.134913, 0.0410350, 0.940134, 0.411034,
            0.484675, 0.169943, 0.325046, 0.269194, 0.552381,
            0.797177, 0.971658, 0.251765, 0.531675, 0.854556,
            0.411237, 0.694380, 0.499562, 0.437242, 0.362451,
            0.343206, 0.524099, 0.158634, 0.728597, 0.198340,
            0.571210, 0.477527, 0.962797, 0.973921, 0.413651,
            0.736380, 0.516366, 0.104283, 0.675993, 0.467053,
            0.230112, 0.866994, 0.469885, 0.964392, 0.541084,
            0.332984, 0.581252, 0.422322, 0.872555, 0.803636,
            0.520998, 0.918942, 0.241564, 0.169263, 0.686649,
            0.708284, 0.707858, 0.00113957, 0.827920, 0.845985,
            0.416961, 0.553842, 0.526549, 0.501051, 0.337514,
            0.700873, 0.152816, 0.762935, 0.650039, 0.483321,
            0.708600, 0.410033, 0.507671, 0.596956, 0.177692,
            0.498112, 0.422037, 0.788333, 0.856578, 0.941245,
            0.432411, 0.356469, 0.341916, 0.0331059, 0.641100,
            0.690452, 0.168667, 0.915178, 0.158406, 0.701508,
            0.841774, 0.434161, 0.153123, 0.420066, 0.0499331,
            0.947241, 0.0768818, 0.410540, 0.843788, 0.0640255,
            0.513463, 0.511104, 0.680434, 0.762480, 0.0563867])
        assert y0.size == 100
        y = smooth(y0, 10)
        assert y.size == 100
        x = np.arange(y0.size, dtype='d')
        sset, outmask = iterfit(x, y, nord=3, maxiter=0, bkspace=10)
        assert sset.npoly == 1
        assert sset.funcname == 'legendre'
        # print(sset)
        # yfit,mask = sset.value(x)
        # print(yfit)
        # pylab.plot(x, y, 'k-', x, yfit, 'r-')
