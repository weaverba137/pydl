# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..spheregroup import spheregroup, spherematch
from .. import PydlutilsException, PydlutilsUserWarning

class TestSpheregroup(object):
    """Test the functions in pydl.pydlutils.spheregroup.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_spheregroup(self, recwarn):
        # np.random.seed(137)
        # Ngroup = 3
        # N = 50
        # spread = 0.05
        linklength = 3.0 # degrees
        # x0 = np.concatenate( ( np.random.normal(loc=1,scale=spread,size=(N,)),
        #                        np.random.normal(loc=-1,scale=spread,size=(N,)),
        #                        np.random.normal(loc=1,scale=spread,size=(N,)),
        #                      )).reshape((N*Ngroup,))
        # y0 = np.concatenate( ( np.random.normal(loc=1,scale=spread,size=(N,)),
        #                        np.random.normal(loc=-1,scale=spread,size=(N,)),
        #                        np.random.normal(loc=-1,scale=spread,size=(N,)),
        #                      )).reshape((N*Ngroup,))
        # z0 = np.concatenate( ( np.random.normal(loc=1,scale=spread,size=(N,)),
        #                        np.random.normal(loc=-1,scale=spread,size=(N,)),
        #                        np.random.normal(loc=0,scale=spread,size=(N,)),
        #                      )).reshape((N*Ngroup,))
        # foo = np.arange(N*Ngroup)
        # np.random.shuffle(foo)
        # x = x0[foo]
        # y = y0[foo]
        # z = z0[foo]
        # r = np.sqrt(x**2 + y**2 + z**2)
        # theta = np.degrees(np.arccos(z/r))
        # phi = np.degrees(np.arctan2(y,x))
        # ra = np.where(phi < 0, phi + 360.0,phi)
        # dec = 90.0 - theta
        # group = spheregroup(ra,dec,linklength)
        #
        # Reproduce IDL results
        #
        ra = np.array([227.48189,44.281521,223.18793,224.92843,41.775477,315.34848,44.242353,316.31759,313.46396,311.91193,43.149445,43.449271,317.18169,224.34194,223.57484,44.701994,228.38290,316.48056,222.82152,224.04855,223.21508,43.901901,49.183680,314.97164,228.35194,316.80579,44.441586,309.97655,313.81441,46.041358,])
        dec = np.array([-34.026420,33.976505,-36.634993,-32.535347,36.290420,0.54110965,35.642976,0.62491526,-0.66827180,-0.083615550,34.695414,34.645375,0.61938086,-35.870473,-35.199912,35.471107,-34.539174,1.5074932,-34.225069,-37.184979,-36.822858,37.650088,35.527187,-1.2731230,-37.131627,-3.1725828,35.996336,1.7463770,-0.95217681,35.159271,])
        expected_ingroup = np.array([0,1,0,0,1,2,1,2,2,2,1,1,2,0,0,1,0,2,0,0,0,1,1,2,0,2,1,2,2,1,])
        expected_multgroup = np.array([10,10,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,])
        expected_firstgroup = np.array([0,1,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,])
        expected_nextgroup = np.array([2,4,3,13,6,7,10,8,9,12,11,15,17,14,16,21,18,23,19,20,24,22,26,25,-1,27,29,28,-1,-1,])
        group = spheregroup(ra,dec,linklength)
        assert (group[0] == expected_ingroup).all()
        assert (group[1] == expected_multgroup).all()
        assert (group[2] == expected_firstgroup).all()
        assert (group[3] == expected_nextgroup).all()
        #
        # Exceptions
        #
        with raises(PydlutilsException):
            group = spheregroup(np.array([137.0]), np.array([55.0]), linklength)
        #
        # warnings
        #
        group = spheregroup(ra, dec, linklength, chunksize=linklength)
        w = recwarn.pop(PydlutilsUserWarning)
        assert "chunksize changed to" in str(w.message)

    def test_spherematch(self):
        i1_should_be = np.array([17, 0, 2, 16, 12, 13, 1, 5, 15, 7,
                                 19, 8, 11, 10, 14, 18, 3, 9, 6, 4])
        i2_should_be = np.array([2, 0, 17, 3, 16, 15, 14, 5, 6, 10,
                                 8, 19, 18, 4, 9, 7, 11, 12, 1, 13])
        np.random.seed(137)
        searchrad = 3.0/3600.0
        n = 20
        ra1 = 360.0*np.random.random((n,))
        dec1 = 90.0 - np.rad2deg(np.arccos(2.0*np.random.random((n,)) - 1.0))
        ra2 = ra1 + np.random.normal(0, 1.0/3600.0)
        dec2 = dec1 + np.random.normal(0, 1.0/3600.0)
        foo = np.arange(n)
        np.random.shuffle(foo)
        i1, i2, d12 = spherematch(ra1, dec1, ra2[foo], dec2[foo], searchrad,
                                  maxmatch=0)
        assert (i1 == i1_should_be).all()
        assert (i2 == i2_should_be).all()
