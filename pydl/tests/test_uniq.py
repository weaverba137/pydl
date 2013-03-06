# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_uniq():
    from ..uniq import uniq
    from numpy import array, sort
    items = array([ 1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1 ])
    u = uniq(sort(items))
    equality = u == array([ 3,  5,  7,  9, 10, 11, 12, 13])
    assert equality.all()
