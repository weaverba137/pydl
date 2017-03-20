# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises, catch_warnings
from astropy.utils.data import get_pkg_data_filename
from ..spheregroup import spheregroup, spherematch
from .. import PydlutilsException, PydlutilsUserWarning


class TestSpheregroup(object):
    """Test the functions in pydl.pydlutils.spheregroup.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_spheregroup(self):
        test_data_file = get_pkg_data_filename('t/spheregroup_data.txt')
        test_data = np.loadtxt(test_data_file, dtype='d', delimiter=',')
        # np.random.seed(137)
        # Ngroup = 3
        # N = 50
        # spread = 0.05
        linklength = 3.0  # degrees
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
        ra = test_data[0, :]
        dec = test_data[1, :]
        expected_ingroup = test_data[2, :].astype(np.int64)
        expected_multgroup = test_data[3, :].astype(np.int64)
        expected_firstgroup = test_data[4, :].astype(np.int64)
        expected_nextgroup = test_data[5, :].astype(np.int64)
        group = spheregroup(ra, dec, linklength)
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
        with catch_warnings(PydlutilsUserWarning) as w:
            group = spheregroup(ra, dec, linklength, chunksize=linklength)
        # w = recwarn.pop(PydlutilsUserWarning)
        assert "chunksize changed to" in str(w[0].message)

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
