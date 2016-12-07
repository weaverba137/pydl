# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import catch_warnings, raises
from .. import PydlutilsUserWarning
from ..rgbcolor import (nw_arcsinh, nw_cut_to_box, nw_float_to_byte,
                        nw_scale_rgb)


class TestRGBColor(object):
    """Test the functions in pydl.pydlutils.rgbcolor.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_nw_arcsinh(self):
        colors = np.random.random((10, 10))
        with raises(ValueError):
            fitted_colors = nw_arcsinh(colors)
        colors = np.random.random((10, 10, 5))
        with raises(ValueError):
            fitted_colors = nw_arcsinh(colors)
        colors = np.random.random((10, 10, 3))
        fitted_colors = nw_arcsinh(colors, nonlinearity=0)
        assert (fitted_colors == colors).all()
        colors = np.ones((2, 2, 3))
        fac = np.arcsinh(9.0)/9.0
        fitted_colors = nw_arcsinh(colors)
        assert np.allclose(fitted_colors, fac)

    def test_nw_cut_to_box(self):
        colors = np.random.random((10, 10))
        with raises(ValueError):
            boxed_colors = nw_cut_to_box(colors)
        colors = np.random.random((10, 10, 5))
        with raises(ValueError):
            boxed_colors = nw_cut_to_box(colors)
        colors = np.random.random((10, 10, 3))
        with raises(ValueError):
            boxed_colors = nw_cut_to_box(colors, origin=(1.0, 1.0))
        boxed_colors = nw_cut_to_box(colors)
        assert np.allclose(boxed_colors, colors)

    def test_nw_float_to_byte(self):
        colors = np.zeros((10, 10, 3), dtype=np.float)
        byte_colors = nw_float_to_byte(colors)
        assert (byte_colors == 0).all()
        colors = np.ones((10, 10, 3), dtype=np.float)
        byte_colors = nw_float_to_byte(colors)
        assert (byte_colors == 255).all()
        with catch_warnings(PydlutilsUserWarning) as w:
            byte_colors = nw_float_to_byte(colors, bits=16)
        assert len(w) > 0

    def test_nw_scale_rgb(self):
        colors = np.random.random((10, 10))
        with raises(ValueError):
            scaled_colors = nw_scale_rgb(colors)
        colors = np.random.random((10, 10, 5))
        with raises(ValueError):
            scaled_colors = nw_scale_rgb(colors)
        colors = np.random.random((10, 10, 3))
        with raises(ValueError):
            scaled_colors = nw_scale_rgb(colors, scales=(1.0, 1.0))
        colors = np.ones((2, 2, 3))
        scaled_colors = nw_scale_rgb(colors, scales=(2.0, 2.0, 2.0))
        assert np.allclose(scaled_colors, 2.0)
