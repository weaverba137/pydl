# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_sdss_calibv():
    from numpy.testing import assert_allclose
    from .. import sdss_calibv
    assert_allclose(0.2650306748466258,sdss_calibv())
