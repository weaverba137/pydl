# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_spheregroup(recwarn):
    import numpy as np
    from astropy.tests.helper import raises
    from .. import spheregroup
    from ... import PydlutilsException, PydlutilsUserWarning
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
        group = spheregroup(np.array([137.0]),np.array([55.0]),linklength)
    #
    # warnings
    #
    group = spheregroup(ra,dec,linklength,chunksize=linklength)
    w = recwarn.pop(PydlutilsUserWarning)
    assert "chunksize changed to" in str(w.message)
