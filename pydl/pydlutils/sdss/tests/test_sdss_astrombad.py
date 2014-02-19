# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from astropy.tests.helper import remote_data
#
@remote_data
def test_sdss_astrombad():
    from .. import sdss_astrombad
    assert sdss_astrombad(77,1,20) == False
    assert sdss_astrombad(77,3,35) == True
    assert sdss_astrombad(77,6,77) == False
