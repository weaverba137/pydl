# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..cooling import read_ds_cooling


class TestCooling(object):
    """Test the functions in pydl.pydlutils.cooling.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_read_ds_cooling(self):
        with raises(ValueError):
            logT, loglambda = read_ds_cooling('m-99.cie')
        logT, logL = read_ds_cooling('m-15.cie')
        assert np.allclose(logT[0:5], np.array([4.0, 4.05, 4.1, 4.15, 4.2]))
        assert np.allclose(logL[0:5], np.array([-26.0, -24.66, -23.52,
                                                -22.62, -22.11]))
        logT = np.array([4.025, 4.125, 4.225, 4.325])
        logT2, logL = read_ds_cooling('m-15.cie', logT)
        assert np.allclose(logT2, logT)
        assert np.allclose(logL, np.array([-25.33, -23.07, -22.04, -22.06]))
