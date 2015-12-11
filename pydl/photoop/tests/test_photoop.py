# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..photoobj import sdss_calibv, unwrap_objid


class TestPhotoobj(object):
    """Test the functions in pydl.photoop.photoobj.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_sdss_calibv(self):
        assert np.allclose(0.2650306748466258, sdss_calibv())

    def test_unwrap_objid(self):
        objid = unwrap_objid(np.array([1237661382772195474]))
        assert objid.skyversion == 2
        assert objid.rerun == 301
        assert objid.run == 3704
        assert objid.camcol == 3
        assert objid.frame == 91
        assert objid.id == 146
        objid = unwrap_objid(np.array(['1237661382772195474']))
        assert objid.skyversion == 2
        assert objid.rerun == 301
        assert objid.run == 3704
        assert objid.camcol == 3
        assert objid.frame == 91
        assert objid.id == 146
        with raises(ValueError):
            objid = unwrap_objid(np.array([3.14159]))
