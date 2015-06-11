# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from .. import sdss_path
from os import environ
from pytest import raises
from .PathTestCase import PathTestCase
#
def TestSdssPath(PathTestCase):
    """Set of tests for sdss_path().
    """
    path_data = {
        'apObj':"/PHOTO_REDUX/301/137/objcs/4",
        'calibMatch':"/PHOTO_REDUX/301/137/nfcalib",
        'calibPhotom':"/PHOTO_REDUX/301/137/nfcalib",
        'calibPhotomGlobal':"/PHOTO_CALIB/301/137/nfcalib",
        'fakeIdR':"/PHOTO_DATA/137/fake_fields/4",
        'fpAtlas':"/PHOTO_REDUX/301/137/objcs/4",
        'fpBIN':"/PHOTO_REDUX/301/137/objcs/4",
        'fpC':"/PHOTO_REDUX/301/137/objcs/4",
        'fpFieldStat':"/PHOTO_REDUX/301/137/objcs/4",
        'fpM':"/PHOTO_REDUX/301/137/objcs/4",
        'fpObjc':"/PHOTO_REDUX/301/137/objcs/4",
        'hoggObj':"/PHOTO_REDUX/301/137/objcs/4",
        'idFF':"/PHOTO_REDUX/301/137/objcs/4",
        'idR':"/PHOTO_DATA/137/fields/4",
        'idRR':"/PHOTO_DATA/137/fields/4",
        'psBB':"/PHOTO_REDUX/301/137/objcs/4",
        'psFF':"/PHOTO_REDUX/301/137/objcs/4",
        'psField':"/PHOTO_REDUX/301/137/objcs/4",
        'tsField':"/PHOTO_REDUX/301/137/calibChunks/4",
        }
    def test_sdss_path_bad_ftype(self):
        #
        # Bad ftype
        #
        with raises(KeyError):
            p = sdss_path('fooBar',137)
    def test_sdss_path(self):
        for ftype in self.path_data:
            assert sdss_path(ftype,137,4,'301') == self.path_data[ftype]
