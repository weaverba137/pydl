# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from os import environ
#
class PathTestCase(object):
    """Set up environment variables for testing sdss_path and sdss_name.
    """
    path_envvars = ('PHOTO_CALIB','PHOTO_DATA','BOSS_PHOTOOBJ','PHOTO_REDUX',
        'PHOTO_RESOLVE','PHOTO_SKY','PHOTO_SWEEP')
    env_cache = {}
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
        'reObjGlobal':"/PHOTO_RESOLVE/301/137/resolve/4",
        'reObjRun':"/PHOTO_REDUX/301/137/resolve/4",
        'reObjTmp':"/PHOTO_RESOLVE/301/137/resolve/4",
        'tsField':"/PHOTO_REDUX/301/137/calibChunks/4",
        }
    name_data = {
        'apObj':"apObj-000137-r4-0042.fit",
        'calibMatch':"calibMatch-000137-4.fits",
        'calibPhotom':"calibPhotom-000137-4.fits",
        'calibPhotomGlobal':"calibPhotomGlobal-000137-4.fits",
        'fakeIdR':"idR-000137-r4-0042.fit",
        'fpAtlas':"fpAtlas-000137-4-0042.fit",
        'fpBIN':"fpBIN-000137-r4-0042.fit",
        'fpC':"fpC-000137-r4-0042.fit",
        'fpFieldStat':"fpFieldStat-000137-4-0042.fit",
        'fpM':"fpM-000137-r4-0042.fit",
        'fpObjc':"fpObjc-000137-4-0042.fit",
        'hoggObj':"hoggObj-000137-4-0042.fits",
        'idFF':"idFF-000137-r4.fit",
        'idR':"idR-000137-r4-0042.fit",
        'idRR':"idRR-000137-r4-0042.fit",
        'psBB':"psBB-000137-r4-0042.fit",
        'psFF':"psFF-000137-r4.fit",
        'psField':"psField-000137-4-0042.fit",
        'reObjGlobal':"reObjGlobal-000137-4-0042.fits",
        'reObjRun':"reObjRun-000137-4-0042.fits",
        'reObjTmp':"reObjTmp-000137-4-0042.fits",
        'tsField':"tsField-000137-4-301-0042.fit",
        }
    def setup(self):
        for e in self.path_envvars:
            try:
                actual_env = environ[e]
            except KeyError:
                actual_env = None
            self.env_cache[e] = actual_env
            environ[e] = '/'+e
        return
    def teardown(self):
        for e in self.env_cache:
            if self.env_cache[e] is None:
                del environ[e]
            else:
                environ[e] = self.env_cache[e]
        return
