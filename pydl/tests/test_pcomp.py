# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


def test_pcomp():
    from os.path import dirname, join
    from ..pcomp import pcomp
    from numpy import array, loadtxt, tile
    try:
        from astropy.tests.compat import assert_allclose
    except ImportError:
        from numpy.testing.utils import assert_allclose
    test_data_file = join(dirname(__file__), 't', 'pcomp_data.txt')
    test_data = loadtxt(test_data_file, dtype='d', delimiter=',')
    pcomp_data = test_data[0:20, :]
    m = 4
    n = 20
    means = tile(pcomp_data.mean(0), 20).reshape(pcomp_data.shape)
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
        assert_allclose(abs(foo.coefficients[:, k]), abs(coefficientsT[:, k]),
                        1e-4)
    eigenvalues = test_data[44, :]
    assert_allclose(foo.eigenvalues, eigenvalues, 1e-4)
    variance = test_data[45, :]
    assert_allclose(foo.variance, variance, 1e-4)
