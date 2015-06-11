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
