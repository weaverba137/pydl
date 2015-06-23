# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from .. import sdss_path
from os import environ
from astropy.tests.helper import raises
from .PathTestCase import PathTestCase
#
class TestSdssPath(PathTestCase):
    """Set of tests for sdss_path().
    """
    def test_sdss_path_bad_ftype(self):
        #
        # Bad ftype
        #
        with raises(KeyError):
            p = sdss_path('fooBar',137)
    def test_sdss_path(self):
        for ftype in self.path_data:
            assert sdss_path(ftype,137,4,'301') == self.path_data[ftype]
