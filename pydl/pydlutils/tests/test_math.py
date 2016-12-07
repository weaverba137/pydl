# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from astropy.utils.data import get_pkg_data_filename
from ..math import computechi2, djs_median, find_contiguous


class TestMath(object):
    """Test the functions in pydl.pydlutils.math.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_computechi2(self):
        x = np.arange(20)
        y = np.array([6.6198438, 1.3491303, 0.41035045, 9.4013375, 4.1103360,
            4.3522868, 4.6338078, 4.7400367, 5.1860726, 5.1082748,
            5.1084596, 5.2990997, 5.5987537, 5.7007504, 5.7855296,
            4.1123709, 6.9437957, 4.9956179, 4.3724215, 3.6245063])
        sqivar = np.array([1.32531, 0.886090, 1.08384, 1.04489, 1.46807,
            1.30800, 0.507725, 1.12840, 0.955025, 1.35925,
            1.10126, 1.45690, 0.575700, 0.949710, 1.23368,
            0.536489, 0.772543, 0.957729, 0.883976, 1.11559])
        templates = np.vstack((np.ones((20,), dtype='d'), x)).transpose()
        chi2 = computechi2(y, sqivar, templates)
        #
        # 20 data points, 2 parameters = 18 degrees of freedom.
        #
        assert chi2.dof == 18
        #
        # The covariance matrix should be symmetric
        #
        assert (chi2.covar.T == chi2.covar).all()
        #
        # The variance vector should be the same as the diagonal of the
        # covariance matrix.
        #
        assert (chi2.var == np.diag(chi2.covar)).all()

    def test_djs_median(self):
        test_data_file = get_pkg_data_filename('t/median_data.txt')
        test_data = np.loadtxt(test_data_file, dtype='d', delimiter=',')
        np.random.seed(424242)
        data = 100.0*np.random.random(100)
        data2 = 100.0*np.random.random((10, 10))
        data3 = 100.0*np.random.random((10, 10, 10))
        data_width_5 = test_data[0, :]
        data_width_5_reflect = test_data[1, :]
        #
        # Degenerate cases that fall back on numpy.median().
        #
        assert np.allclose(np.median(data), djs_median(data))
        assert np.allclose(np.median(data2, axis=0),
                            djs_median(data2, dimension=0))
        assert np.allclose(np.median(data2, axis=1),
                            djs_median(data2, dimension=1))
        #
        # Test widths.
        #
        assert np.allclose(data, djs_median(data, width=1))
        assert np.allclose(data_width_5, djs_median(data, width=5))
        # assert np.allclose(data_width_5_reflect,
        #                 djs_median(data, width=5, boundary='reflect'))
        #
        # Exceptions
        #
        with raises(ValueError):
            foo = djs_median(data, width=5, boundary='nearest')
        with raises(ValueError):
            foo = djs_median(data, width=5, boundary='wrap')
        with raises(ValueError):
            foo = djs_median(data, width=5, boundary='foobar')
        with raises(ValueError):
            foo = djs_median(data2, width=5, dimension=1)
        with raises(ValueError):
            foo = djs_median(data3, width=5)
        with raises(ValueError):
            foo = djs_median(data3, width=5, boundary='reflect')

    def test_find_contiguous(self):
        assert (find_contiguous(np.array([0, 1, 1, 1, 0, 1, 1, 0, 1])) ==
                [1, 2, 3])
