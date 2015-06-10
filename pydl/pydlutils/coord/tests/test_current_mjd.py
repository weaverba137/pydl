# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_current_mjd():
    from .. import current_mjd
    assert current_mjd() > 50000.0
