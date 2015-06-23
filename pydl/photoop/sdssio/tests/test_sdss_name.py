# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from .. import sdss_name
from os import environ
from os.path import join
from astropy.tests.helper import raises
from .PathTestCase import PathTestCase
#
class TestSdssName(PathTestCase):
    """Set of tests for sdss_path().
    """
    def test_sdss_name_bad_ftype(self):
        #
        # Bad ftype
        #
        with raises(KeyError):
            p = sdss_name('fooBar',137,4,42)
    def test_sdss_name(self):
        for ftype in self.name_data:
            assert sdss_name(ftype,137,4,42,'301','r',no_path=True) == self.name_data[ftype]
    def test_sdss_name_with_path(self):
        for ftype in self.name_data:
            assert sdss_name(ftype,137,4,42,'301','r') == join(self.path_data[ftype],self.name_data[ftype])
    def test_sdss_name_reObj(self):
        assert sdss_name('reObj',137,4,42,'301','r',no_path=True) == self.name_data['reObjGlobal']
        resolve = environ['PHOTO_RESOLVE']
        del environ['PHOTO_RESOLVE']
        assert sdss_name('reObj',137,4,42,'301','r',no_path=True) == self.name_data['reObjRun']
        environ['PHOTO_RESOLVE'] = resolve
