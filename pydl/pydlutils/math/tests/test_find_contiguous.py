# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_find_contiguous():
    from numpy import array
    from .. import find_contiguous
    assert find_contiguous(array([0,1,1,1,0,1,1,0,1])) == [1,2,3]
