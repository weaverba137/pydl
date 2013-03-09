# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_find_contiguous():
    from .. import sdss_objid
    assert sdss_objid(3704,3,91,146) == 1237661382772195474
