# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_sdss_calib():
    from .. import sdss_calib
    foo = sdss_calib(94,6,101)
    assert foo['NMGYPERCOUNT'] == 1.0
