# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_unwrap_objid():
    from numpy import array
    from astropy.tests.helper import raises
    from .. import unwrap_objid
    objid = unwrap_objid(array([1237661382772195474]))
    assert objid.skyversion == 2
    assert objid.rerun == 301
    assert objid.run == 3704
    assert objid.camcol == 3
    assert objid.frame == 91
    assert objid.id == 146
    objid = unwrap_objid(array(['1237661382772195474']))
    assert objid.skyversion == 2
    assert objid.rerun == 301
    assert objid.run == 3704
    assert objid.camcol == 3
    assert objid.frame == 91
    assert objid.id == 146
    with raises(ValueError):
        objid = unwrap_objid(array([3.14159]))
