# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_sdss_objid():
    from numpy import array
    from astropy.tests.helper import raises
    from .. import sdss_objid
    assert sdss_objid(3704,3,91,146) == 1237661382772195474
    run = array([3704,1000])
    camcol = array([3,6])
    field = array([91,77])
    obj = array([146,123])
    assert (array([1237661382772195474,1237649770790322299]) == sdss_objid(run,camcol,field,obj)).all()
    #
    # Exceptions
    #
    with raises(ValueError):
        objid = sdss_objid(run,3,91,146)
    with raises(ValueError):
        objid = sdss_objid(3704,camcol,91,146)
    with raises(ValueError):
        objid = sdss_objid(3704,3,field,146)
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,obj)
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,146,rerun=array([137,301]))
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,146,skyversion=array([2,3]))
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,146,skyversion=-2)
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,146,skyversion=16)
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,146,rerun=-2)
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,146,rerun=2**11)
    with raises(ValueError):
        objid = sdss_objid(-2,3,91,146)
    with raises(ValueError):
        objid = sdss_objid(2**16,3,91,146)
    with raises(ValueError):
        objid = sdss_objid(3704,0,91,146)
    with raises(ValueError):
        objid = sdss_objid(3704,7,91,146)
    with raises(ValueError):
        objid = sdss_objid(3704,3,-2,146)
    with raises(ValueError):
        objid = sdss_objid(3704,3,2**12,146)
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,-2)
    with raises(ValueError):
        objid = sdss_objid(3704,3,91,2**16)
