# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test the top-level pydl functions.
"""
import pytest
import numpy as np
try:
    from astropy.tests.compat import assert_allclose
except ImportError:
    from numpy.testing import assert_allclose
from astropy.tests.helper import raises
from astropy.utils.data import get_pkg_data_filename
from ..file_lines import file_lines
from ..median import median
from ..pcomp import pcomp
from ..rebin import rebin
from ..smooth import smooth
from ..uniq import uniq


def test_file_lines():
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


def test_median():
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


def test_pcomp():
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


def test_rebin():
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


@pytest.mark.xfail
def test_rebin_int():
    """Test rebin on integers.  Comparing to IDL code similar to this::

        IDL> seed = 100
        IDL> array = FIX(RANDOMN(seed, 10, 20) * 100, TYPE=1) ; UINT8
        IDL> array_rebin = REBIN(array, 5, 10)
    """
    array = np.array([[188, 186,  25, 212,  34,  98,   3, 235, 155, 148],
                      [107, 166,   4,  41, 101, 190,  39, 154, 153, 239],
                      [135, 181,  92, 161, 213, 136,  35,  61,  80, 164],
                      [123, 248,   8, 157,  96, 118,  99,   1, 109, 246],
                      [226,  71, 183,  27,  46,  99,   8, 239,  66,  25],
                      [ 27, 219,  37, 130,   5,  81,  65, 250,  96,  14],
                      [ 71, 157, 156, 136,  47, 225, 247, 191,  49,  12],
                      [231, 133,   9,  38, 243,   2, 235, 145,  23,  22],
                      [146,  38,  49,  89,  42,  57, 220, 214, 135,  47],
                      [101, 116, 122, 209, 141,  37, 158, 224, 245,  82],
                      [ 15,  47,  51, 250, 207, 193, 209, 228, 110,   1],
                      [ 59, 232, 216, 224,  24, 118, 190,  10, 107,  27],
                      [ 84, 193, 112, 206, 113, 171, 138, 117, 244,  20],
                      [  5,  31, 128, 214, 200, 119,  59,  27,  57,  10],
                      [226,  71, 177,  85,   0,  68,  54, 207, 141, 250],
                      [ 52, 119, 121, 177, 165,  99,  68,  29, 137, 200],
                      [172,  91, 181, 187,  87, 250,  45, 154,  58,  83],
                      [ 56, 175, 189,  35, 203, 223, 243, 187, 252,  97],
                      [186, 172, 207, 128,  61, 231,  89,  57, 131, 222],
                      [206,  96,  29,  60,   3,   8, 221,  55,  60,  17]], dtype=np.uint8)
    array_rebin = np.array([[161,  70, 105, 107, 173],
                            [171, 104, 140,  49, 149],
                            [135,  94,  57, 140,  50],
                            [148,  84, 129, 204,  26],
                            [100, 117,  69, 204, 127],
                            [ 88, 185, 135, 159,  61],
                            [ 78, 165, 150,  85,  82],
                            [116, 140,  83,  89, 181],
                            [123, 148, 190, 157, 122],
                            [165, 105,  75, 105, 107]], dtype=np.uint8)

    ar = rebin(array, (10, 5))
    assert (array_rebin == ar).all()


def test_rebin_int_sample():
    """Similar to test_rebin_int(), but using the sample option.
    """
    array = np.array([[188, 186,  25, 212,  34,  98,   3, 235, 155, 148],
                      [107, 166,   4,  41, 101, 190,  39, 154, 153, 239],
                      [135, 181,  92, 161, 213, 136,  35,  61,  80, 164],
                      [123, 248,   8, 157,  96, 118,  99,   1, 109, 246],
                      [226,  71, 183,  27,  46,  99,   8, 239,  66,  25],
                      [ 27, 219,  37, 130,   5,  81,  65, 250,  96,  14],
                      [ 71, 157, 156, 136,  47, 225, 247, 191,  49,  12],
                      [231, 133,   9,  38, 243,   2, 235, 145,  23,  22],
                      [146,  38,  49,  89,  42,  57, 220, 214, 135,  47],
                      [101, 116, 122, 209, 141,  37, 158, 224, 245,  82],
                      [ 15,  47,  51, 250, 207, 193, 209, 228, 110,   1],
                      [ 59, 232, 216, 224,  24, 118, 190,  10, 107,  27],
                      [ 84, 193, 112, 206, 113, 171, 138, 117, 244,  20],
                      [  5,  31, 128, 214, 200, 119,  59,  27,  57,  10],
                      [226,  71, 177,  85,   0,  68,  54, 207, 141, 250],
                      [ 52, 119, 121, 177, 165,  99,  68,  29, 137, 200],
                      [172,  91, 181, 187,  87, 250,  45, 154,  58,  83],
                      [ 56, 175, 189,  35, 203, 223, 243, 187, 252,  97],
                      [186, 172, 207, 128,  61, 231,  89,  57, 131, 222],
                      [206,  96,  29,  60,   3,   8, 221,  55,  60,  17]], dtype=np.uint8)
    array_sample = np.array([[188,  25,  34,   3, 155],
                             [135,  92, 213,  35,  80],
                             [226, 183,  46,   8,  66],
                             [ 71, 156,  47, 247,  49],
                             [146,  49,  42, 220, 135],
                             [ 15,  51, 207, 209, 110],
                             [ 84, 112, 113, 138, 244],
                             [226, 177,   0,  54, 141],
                             [172, 181,  87,  45,  58],
                             [186, 207,  61,  89, 131]], dtype=np.uint8)
    ars = rebin(array, (10, 5), sample=True)
    assert (array_sample == ars).all()


def test_smooth():
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


def test_uniq():
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
