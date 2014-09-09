# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_filternum():
    from .. import filternum
    assert filternum('u') == 0
    assert filternum('g') == 1
    assert filternum('r') == 2
    assert filternum('i') == 3
    assert filternum('z') == 4
    #
    # Test default return value
    #
    fn = filternum()
    for k in range(5):
        assert fn[k] == k
