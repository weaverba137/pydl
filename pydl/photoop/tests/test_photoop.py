# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test the functions in pydl.photoop.photoobj.
"""
import pytest
import numpy as np
from ..photoobj import sdss_calibv, unwrap_objid


def test_sdss_calibv():
    assert np.allclose(0.2650306748466258, sdss_calibv().value)


def test_unwrap_objid():
    objid = unwrap_objid(np.array([1237661382772195474]))
    assert objid.skyversion == 2
    assert objid.rerun == 301
    assert objid.run == 3704
    assert objid.camcol == 3
    assert objid.firstfield == 0
    assert objid.frame == 91
    assert objid.id == 146
    objid = unwrap_objid(np.array(['1237661382772195474']))
    assert objid.skyversion == 2
    assert objid.rerun == 301
    assert objid.run == 3704
    assert objid.camcol == 3
    assert objid.firstfield == 0
    assert objid.frame == 91
    assert objid.id == 146
    objid = unwrap_objid(np.array([587722984180548043]))
    assert objid.skyversion == 1
    assert objid.rerun == 40
    assert objid.run == 752
    assert objid.camcol == 5
    assert objid.firstfield == 1
    assert objid.frame == 618
    assert objid.id == 459
    with pytest.raises(ValueError):
        objid = unwrap_objid(np.array([3.14159]))
