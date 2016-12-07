# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
import glob
try:
    from astropy.tests.compat import assert_allclose
except ImportError:
    from numpy.testing.utils import assert_allclose
from astropy.tests.helper import raises
from astropy.utils.data import get_pkg_data_filename
from os.path import basename, dirname, join
from ..file_lines import file_lines
from ..median import median
from ..pcomp import pcomp
from ..rebin import rebin
from ..smooth import smooth
from ..uniq import uniq


class TestPydl(object):
    """Test the top-level pydl functions.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_file_lines(self):
        #
        # Find the test files
        #
        line_numbers = (1, 42, 137)
        plainfiles = [get_pkg_data_filename('t/this-file-contains-{0:d}-lines.txt'.format(l)) for l in line_numbers]
        gzfiles = [get_pkg_data_filename('t/this-file-contains-{0:d}-lines.txt.gz'.format(l)) for l in line_numbers]
        for i, p in enumerate(plainfiles):
            n = file_lines(p)
            assert n == line_numbers[i]
        for i, p in enumerate(gzfiles):
            n = file_lines(p, compress=True)
            assert n == line_numbers[i]
        #
        # Test list passing
        #
        n = file_lines(plainfiles)
        assert tuple(n) == line_numbers
        n = file_lines(gzfiles, compress=True)
        assert tuple(n) == line_numbers
        #
        # Make sure empty files work
        #
        n = file_lines(get_pkg_data_filename('t/this-file-is-empty.txt'))
        assert n == 0

    def test_median(self):
        odd_data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                            dtype=np.float32)
        even_data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                             dtype=np.float32)
        assert median(odd_data) == 7
        assert median(odd_data, even=True) == 7
        assert median(even_data) == 7
        assert median(even_data, even=True) == 6.5
        assert (median(odd_data, 3) == odd_data).all()
        with raises(ValueError):
            foo = median(np.ones((9, 9, 9)), 3)
        odd_data2 = np.vstack((odd_data, odd_data, odd_data, odd_data, odd_data))
        assert (median(odd_data2, 3) == odd_data2).all()
        assert (median(odd_data2, axis=0) == odd_data).all()
        assert (median(odd_data2, axis=1) ==
                7*np.ones((odd_data2.shape[0],), dtype=odd_data2.dtype)).all()

    def test_pcomp(self):
        test_data_file = get_pkg_data_filename('t/pcomp_data.txt')
        test_data = np.loadtxt(test_data_file, dtype='d', delimiter=',')
        with raises(ValueError):
            foo = pcomp(np.arange(10))
        pcomp_data = test_data[0:20, :]
        m = 4
        n = 20
        means = np.tile(pcomp_data.mean(0), n).reshape(pcomp_data.shape)
        newarray = pcomp_data - means
        foo = pcomp(newarray, covariance=True)
        #
        # This array is obtained from the IDL version of PCOMP.
        # It is only accurate up to an overall sign on each column.
        #
        derived = test_data[20:40, :]
        for k in range(m):
            assert_allclose(abs(foo.derived[:, k]), abs(derived[:, k]), 1e-4)
        coefficients = test_data[40:44, :]
        coefficientsT = coefficients.T
        for k in range(m):
            assert_allclose(abs(foo.coefficients[:, k]),
                            abs(coefficientsT[:, k]),
                            1e-4)
        eigenvalues = test_data[44, :]
        assert_allclose(foo.eigenvalues, eigenvalues, 1e-4)
        variance = test_data[45, :]
        assert_allclose(foo.variance, variance, 1e-4)
        #
        # Test the standardization.
        #
        foo = pcomp(pcomp_data, standardize=True, covariance=True)
        # for k in range(m):
        #     assert_allclose(abs(foo.derived[:, k]), abs(derived[:, k]), 1e-4)
        # for k in range(m):
        #     assert_allclose(abs(foo.coefficients[:, k]),
        #                     abs(coefficientsT[:, k]),
        #                     1e-4)
        eigenvalues = test_data[46, :]
        assert_allclose(foo.eigenvalues, eigenvalues, 1e-4)
        variance = test_data[47, :]
        assert_allclose(foo.variance, variance, 1e-4)
        # assert_allclose(foo.derived[0, :], np.array([-1.64153312,
        #                                              -9.12322038,
        #                                              1.41790708,
        #                                              -8.29359322]))
        #
        # Make sure correlation is working at least.
        #
        foo = pcomp(pcomp_data, standardize=True)
        assert_allclose(foo.eigenvalues, np.array([2.84968632e+00,
                                                   1.00127640e+00,
                                                   1.48380121e-01,
                                                   6.57156222e-04]))
        assert_allclose(foo.variance, np.array([7.12421581e-01,
                                                2.50319100e-01,
                                                3.70950302e-02,
                                                1.64289056e-04]))

    def test_rebin(self):
        x = np.arange(40)
        with raises(ValueError):
            r = rebin(x, d=(10, 10))
        with raises(ValueError):
            r = rebin(x, d=(70,))
        with raises(ValueError):
            r = rebin(x, d=(30,))
        x = np.array([[1.0, 2.0], [2.0, 3.0]])
        rexpect = np.array([[1.0, 2.0], [1.5, 2.5], [2.0, 3.0], [2.0, 3.0]])
        r = rebin(x, d=(4, 2))
        assert np.allclose(r, rexpect)
        rexpect = np.array([[1.0, 1.5, 2.0, 2.0], [2.0, 2.5, 3.0, 3.0]])
        r = rebin(x, d=(2, 4))
        assert np.allclose(r, rexpect)
        rexpect = np.array([[1.0, 2.0], [1.0, 2.0], [2.0, 3.0], [2.0, 3.0]])
        r = rebin(x, d=(4, 2), sample=True)
        assert np.allclose(r, rexpect)
        rexpect = np.array([[1.0, 1.0, 2.0, 2.0], [2.0, 2.0, 3.0, 3.0]])
        r = rebin(x, d=(2, 4), sample=True)
        assert np.allclose(r, rexpect)
        x = np.arange(10)
        rexpect = np.array([0.0, 2.0, 4.0, 6.0, 8.0])
        r = rebin(x, d=(5,), sample=True)
        assert np.allclose(r, rexpect)
        x = np.array([[1.0, 2.0, 3.0, 4.0],
                      [2.0, 3.0, 4.0, 5.0],
                      [3.0, 4.0, 5.0, 6.0],
                      [4.0, 5.0, 6.0, 7.0]])
        rexpect = np.array([[2.0, 4.0], [4.0, 6.0]])
        r = rebin(x, d=(2, 2))
        assert np.allclose(r, rexpect)
        rexpect = np.array([[1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 4.5],
                            [3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 6.5]])
        r = rebin(x, d=(2, 8))
        assert np.allclose(r, rexpect)

    def test_smooth(self):
        test_data_file = get_pkg_data_filename('t/smooth_data.txt')
        noise = np.loadtxt(test_data_file, dtype='d')
        #
        # Test smooth function
        #
        x = 8.0*np.arange(100)/100.0 - 4.0
        y = np.sin(x) + 0.1*noise
        s = smooth(y, 5)
        assert s.shape == (100,)
        s_edge = smooth(y, 5, True)
        assert s_edge.shape == (100,)
        s_w = smooth(y, 1)
        assert (s_w == y).all()

    def test_uniq(self):
        items = np.array([1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1])
        items_sorted = np.sort(items)
        items_argsorted = np.argsort(items)
        #
        # Test pre-sorted array.
        #
        u1 = uniq(items_sorted)
        assert (u1 == np.array([3, 5, 7, 9, 10, 11, 12, 13])).all()
        #
        # Test arg-sorted array.
        #
        u2 = uniq(items, items_argsorted)
        assert (u2 == np.array([13, 9, 8, 10, 5, 7, 11, 12])).all()
        assert (items_sorted[u1] == items[u2]).all()
        #
        # Test degenerate case of all identical items.
        #
        identical_items = np.ones((10,), dtype=items.dtype)
        u = uniq(identical_items)
        assert (u == np.array([9])).all()
        u = uniq(identical_items, np.arange(10, dtype=items.dtype))
        assert (u == np.array([9])).all()
