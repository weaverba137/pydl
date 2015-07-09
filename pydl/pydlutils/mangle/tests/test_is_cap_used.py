# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_is_cap_used():
    from .. import is_cap_used
    assert is_cap_used(1<<2, 2)
    assert not is_cap_used(1<<2,1)
