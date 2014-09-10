# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
#
def test_uniq():
    from ..uniq import uniq
    from numpy import array, sort, argsort, ones, arange
    items = array([ 1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1 ])
    items_sorted = sort(items)
    items_argsorted = argsort(items)
    #
    # Test pre-sorted array.
    #
    u1 = uniq(items_sorted)
    equality = u1 == array([ 3,  5,  7,  9, 10, 11, 12, 13])
    assert equality.all()
    #
    # Test arg-sorted array.
    #
    u2 = uniq(items,items_argsorted)
    equality = u2 == array([ 13,  9,  8,  10, 5, 7, 11, 12])
    assert equality.all()
    assert (items_sorted[u1] == items[u2]).all()
    #
    # Test degenerate case of all identical items.
    #
    identical_items = ones((10,),dtype=items.dtype)
    u = uniq(identical_items)
    equality = u == array([9])
    assert equality.all()
    u = uniq(identical_items,arange(10,dtype=items.dtype))
