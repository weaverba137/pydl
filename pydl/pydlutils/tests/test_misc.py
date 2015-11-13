# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..misc import djs_laxisgen, djs_laxisnum, hogg_iau_name, struct_print


class TestMisc(object):
    """Test the functions in pydl.pydlutils.misc.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_djs_laxisgen(self):
        #
        # 1d
        #
        assert (np.arange(4, dtype='i4') == djs_laxisgen((4,))).all()
        #
        # 2d
        #
        l = np.array([[0, 0, 0, 0],
                        [1, 1, 1, 1],
                        [2, 2, 2, 2],
                        [3, 3, 3, 3]],
                        dtype='i4')
        assert (l == djs_laxisgen((4, 4))).all()
        assert (l.T == djs_laxisgen((4, 4), iaxis=1)).all()
        with raises(ValueError):
            foo = djs_laxisgen((4, 4), iaxis=2)
        #
        # 3d
        #
        l = np.zeros((4, 4, 4), dtype='i4')
        l[1, :, :] = 1
        l[2, :, :] = 2
        l[3, :, :] = 3
        assert (l == djs_laxisgen((4, 4, 4))).all()
        assert (l.swapaxes(0, 1) == djs_laxisgen((4, 4, 4), iaxis=1)).all()
        assert (l.swapaxes(0, 2) == djs_laxisgen((4, 4, 4), iaxis=2)).all()
        with raises(ValueError):
            foo = djs_laxisgen((4, 4, 4), iaxis=3)
        #
        # More d
        #
        with raises(ValueError):
            foo = djs_laxisgen((4, 4, 4, 4))

    def test_djs_laxisnum(self):
        #
        # 1d
        #
        assert (np.zeros((4,), dtype='i4') == djs_laxisnum((4,))).all()
        #
        # 2d
        #
        l = np.array([[0, 0, 0, 0],
                        [1, 1, 1, 1],
                        [2, 2, 2, 2],
                        [3, 3, 3, 3]],
                        dtype='i4')
        assert (l == djs_laxisnum((4, 4))).all()
        assert (l.T == djs_laxisnum((4, 4), iaxis=1)).all()
        with raises(ValueError):
            foo = djs_laxisnum((4, 4), iaxis=2)
        #
        # 3d
        #
        l = np.zeros((4, 4, 4), dtype='i4')
        l[1, :, :] = 1
        l[2, :, :] = 2
        l[3, :, :] = 3
        assert (l == djs_laxisnum((4, 4, 4))).all()
        assert (l.swapaxes(0, 1) == djs_laxisnum((4, 4, 4), iaxis=1)).all()
        assert (l.swapaxes(0, 2) == djs_laxisnum((4, 4, 4), iaxis=2)).all()
        with raises(ValueError):
            foo = djs_laxisnum((4, 4, 4), iaxis=3)
        #
        # More d
        #
        with raises(ValueError):
            foo = djs_laxisnum((4, 4, 4, 4))

    def test_hogg_iau_name(self):
        assert (hogg_iau_name(354.120375, -0.544777778) ==
                'SDSS J233628.89-003241.2')
        assert (hogg_iau_name(354.120375, -0.544777778, prefix='2MASS') ==
                '2MASS J233628.89-003241.2')
        assert (hogg_iau_name(354.120375, -0.544777778, prefix='') ==
                'J233628.89-003241.2')
        assert (hogg_iau_name(354.120375, -0.544777778, precision=0) ==
                'SDSS J233628.8-003241')
        assert (hogg_iau_name(354.120375, -0.544777778, precision=2) ==
                'SDSS J233628.890-003241.20')
        ra = np.array([354.120375, 7.89439, 36.31915, 110.44730])
        dec = np.array([-0.544777778, -0.35157, 0.47505, 39.35352])
        names = hogg_iau_name(ra, dec)
        assert tuple(names) == ('SDSS J233628.89-003241.2',
                                'SDSS J003134.65-002105.6',
                                'SDSS J022516.59+002830.1',
                                'SDSS J072147.35+392112.6')

    def test_struct_print(self):
        n = 20
        slist = np.zeros(n, dtype=[('PLATE', 'i4'), ('MJD', 'i4'),
                                    ('FIBERID', 'i4'), ('RA', 'f8'),
                                    ('DEC', 'f8'), ('MATCHRAD', 'f4'),
                                    ('RERUN', 'S8')])
        slist['RERUN'][0:int(n/2)] = 'v1_2_3'
        slist['RERUN'][int(n/2):n] = 'v1_3_34'
        slist['PLATE'] = np.random.random_integers(10, 10000, (n,))
        slist['MJD'] = np.random.random_integers(51000, 56000, (n,))
        slist['FIBERID'] = np.random.random_integers(1, 1000, (n,))
        slist['RA'] = 360.0*np.random.random((n,))
        slist['DEC'] = (90.0 -
                        np.rad2deg(np.arccos(
                        2.0*np.random.random((n,)) - 1.0)))
        slist['MATCHRAD'] = np.random.random((n,))
        # print(slist)
        # lines, css = struct_print(slist, debug=True)
        # print(lines)
        # print(css)
        # lines, css = struct_print(slist, debug=True, html=True)
        # print(lines)
        # print(css)
        assert slist.size == n
