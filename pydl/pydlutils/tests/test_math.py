# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
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
        np.random.seed(424242)
        data = 100.0*np.random.random(100)
        data2 = 100.0*np.random.random((10, 10))
        data3 = 100.0*np.random.random((10, 10, 10))
        data_width_5 = np.array([
            49.68311576, 74.83671757, 49.68311576, 42.12573137,
            32.11576752, 27.22092569, 15.08905903, 15.08905903,
            27.22092569, 20.07809411, 41.55978953, 41.55978953,
            35.01184343, 35.01184343, 38.7967727,  38.7967727,
            38.7967727,  38.7967727,  38.7967727,  28.85780425,
            28.85780425, 30.23801035, 30.23801035, 27.0687078,
            30.23801035, 64.77058052, 29.42407045, 38.51749618,
            62.48793217, 38.51749618, 29.42407045, 38.51749618,
            18.46919414, 18.46919414, 42.72395431, 54.11584807,
            54.11584807, 75.67593452, 75.67593452, 72.4380356,
            61.68761955, 61.68761955, 45.36796557, 45.36796557,
            61.68761955, 76.57621138, 76.57621138, 76.57621138,
            23.28589488, 23.28589488, 13.82755909, 12.10607597,
            13.82755909, 27.51891089, 44.21266068, 44.21266068,
            44.21266068, 47.18083025, 47.18083025, 47.18083025,
            47.18083025, 47.18083025, 35.50809933, 35.50809933,
            25.52222293, 25.52222293, 67.8568479,  88.54983822,
            67.8568479,  93.40053148, 93.40053148, 64.12514945,
            47.82074715, 47.82074715, 47.82074715, 34.82234113,
            34.82234113, 52.91092248, 78.51075522, 92.16442338,
            92.16442338, 92.16442338, 72.07989558, 72.07989558,
            68.72579501, 72.07989558, 72.07989558, 70.43134908,
            34.55273356, 62.09010468, 62.09010468, 70.43134908,
            68.89705132, 68.89705132, 68.89705132, 66.30426084,
            55.92748086, 55.92748086, 55.92748086, 65.11387444])
        data_width_5_reflect = np.array([
            49.68311576, 49.68311576, 49.68311576, 42.12573137,
            32.11576752, 27.22092569, 15.08905903, 15.08905903,
            27.22092569, 20.07809411, 41.55978953, 41.55978953,
            35.01184343, 35.01184343, 38.7967727,  38.7967727,
            38.7967727,  38.7967727,  38.7967727,  28.85780425,
            28.85780425, 30.23801035, 30.23801035, 27.0687078,
            30.23801035, 64.77058052, 29.42407045, 38.51749618,
            62.48793217, 38.51749618, 29.42407045, 38.51749618,
            18.46919414, 18.46919414, 42.72395431, 54.11584807,
            54.11584807, 75.67593452, 75.67593452, 72.4380356,
            61.68761955, 61.68761955, 45.36796557, 45.36796557,
            61.68761955, 76.57621138, 76.57621138, 76.57621138,
            23.28589488, 23.28589488, 13.82755909, 12.10607597,
            13.82755909, 27.51891089, 44.21266068, 44.21266068,
            44.21266068, 47.18083025, 47.18083025, 47.18083025,
            47.18083025, 47.18083025, 35.50809933, 35.50809933,
            25.52222293, 25.52222293, 67.8568479,  88.54983822,
            67.8568479,  93.40053148, 93.40053148, 64.12514945,
            47.82074715, 47.82074715, 47.82074715, 34.82234113,
            34.82234113, 52.91092248, 78.51075522, 92.16442338,
            92.16442338, 92.16442338, 72.07989558, 72.07989558,
            68.72579501, 72.07989558, 72.07989558, 70.43134908,
            34.55273356, 62.09010468, 62.09010468, 70.43134908,
            68.89705132, 68.89705132, 68.89705132, 55.92748086,
            65.11387444, 55.92748086, 55.92748086, 55.92748086])
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
