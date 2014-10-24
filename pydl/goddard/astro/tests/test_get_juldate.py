# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_get_juldate():
    from .. import get_juldate
    now = get_juldate()
    assert now > 2400000.5
    assert get_juldate(0) == 40587.0 + 2400000.5
    assert get_juldate(86400) == 40588.0 + 2400000.5
