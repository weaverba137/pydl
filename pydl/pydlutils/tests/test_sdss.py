# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
import pydl.pydlutils.sdss
from astropy.tests.helper import remote_data, raises
from astropy.utils.data import get_pkg_data_filename
from ..sdss import (default_skyversion, sdss_flagexist, sdss_flagname,
                    sdss_flagval, set_maskbits, sdss_astrombad,
                    sdss_objid, sdss_specobjid, unwrap_specobjid)


class TestSDSS(object):
    """Test the functions in pydl.pydlutils.sdss.
    """

    def setup(self):
        pydl.pydlutils.sdss.maskbits = set_maskbits(
            maskbits_file=get_pkg_data_filename('t/testMaskbits.par'))
        self.opbadfields = np.array([
            (77, 'astrom', 30, 73, 'Large astrometric offset at field 39... 72'),
            (85, 'astrom', 8, 28, 'Large astrometric offset at field 11... 27'),
            (85, 'rotator', 242, 253, 'Large rotator offset at field 251 252'),
            (209, 'astrom', 8, 116, 'Tel. offsets before r-band field 115 -DJS'),
            (209, 'astrom', 137, 175, 'Tel. offsets after r-band field 145 -DJS'),
            (250, 'astrom', 456, 468, 'Large astrometric offset -Manual'),
            (251, 'astrom', 147, 159, 'Large astrometric offset at field 156 158')],
            dtype=[('run', '<i4'), ('problem', 'S8'), ('firstfield', '<i4'), ('lastfield', '<i4'), ('comments', 'S47')])
        pydl.pydlutils.sdss.opbadfields = self.opbadfields
        return

    def teardown(self):
        pass

    def test_default_skyversion(self):
        assert default_skyversion() == 2

    def test_sdss_flagexist(self):
        assert sdss_flagexist('TARGET', 'ROSAT_A')
        assert sdss_flagexist('TARGET',
                            ['rosat_a', 'rosat_b', 'rosat_c', 'rosat_d'])
        l, f = sdss_flagexist('target', 'FOOBAR', flagexist=True)
        assert not l
        assert f
        l, which = sdss_flagexist('TARGET', ['rosat_a', 'rosat_b', 'rosat_c',
                                'rosat_d', 'foobar'], whichexist=True)
        assert not l
        assert tuple(which) == (True, True, True, True, False)
        l, f, which = sdss_flagexist('TARGET', ['rosat_a', 'rosat_b',
                                    'rosat_c', 'rosat_d', 'foobar'],
                                    flagexist=True, whichexist=True)
        assert not l
        assert f
        assert tuple(which) == (True, True, True, True, False)

    def test_sdss_flagname(self):
        names = sdss_flagname('ANCILLARY_TARGET1', 2310346608843161600)
        assert tuple(names) == ('BRIGHTGAL', 'BLAZGX', 'ELG')
        names = sdss_flagname('ANCILLARY_TARGET1', 2310346608843161600,
                            concat=True)
        assert names == 'BRIGHTGAL BLAZGX ELG'
        with raises(KeyError):
            names = sdss_flagname('ABADMASK', 123456789)

    def test_sdss_flagval(self):
        val = sdss_flagval('TARGET', 'ROSAT_A')
        assert val == 2**9
        val = sdss_flagval('ANCILLARY_TARGET1', ['BLAZGX', 'ELG', 'BRIGHTGAL'])
        assert val == 2310346608843161600
        with raises(KeyError):
            val = sdss_flagval('TARGET', 'ROSAT_Q')
        with raises(KeyError):
            val = sdss_flagval('ABADMASK', "ABADFLAG")

    def test_set_maskbits(self):
        maskbits = pydl.pydlutils.sdss.maskbits
        assert (set(maskbits.keys()) ==
            set(['TARGET', 'BOSS_TARGET1', 'PRIMTARGET', 'ANCILLARY_TARGET1',
                'ZWARNING', 'TTARGET', 'SECTARGET', 'LEGACY_TARGET2',
                'LEGACY_TARGET1', 'SPECIAL_TARGET2', 'FLUXMATCH_STATUS']))
        assert (set(maskbits['TARGET'].keys()) ==
            set(['QSO_FIRST_SKIRT', 'QSO_CAP', 'GALAXY_RED', 'STAR_CARBON',
                'STAR_WHITE_DWARF', 'GALAXY_RED_II', 'GALAXY_BIG',
                'GALAXY_BRIGHT_CORE', 'SERENDIP_MANUAL', 'STAR_SUB_DWARF',
                'QSO_FIRST_CAP', 'QSO_SKIRT', 'STAR_PN', 'STAR_BHB', 'QSO_HIZ',
                'STAR_BROWN_DWARF', 'SERENDIP_FIRST', 'SOUTHERN_SURVEY',
                'STAR_RED_DWARF', 'STAR_CATY_VAR', 'QSO_REJECT', 'GALAXY',
                'SERENDIP_RED', 'SERENDIP_DISTANT', 'QSO_MAG_OUTLIER',
                'ROSAT_A', 'ROSAT_C', 'ROSAT_B', 'ROSAT_E', 'ROSAT_D',
                'SERENDIP_BLUE']))

    def test_sdss_astrombad(self):
        assert not sdss_astrombad(77, 1, 20)
        assert sdss_astrombad(77, 3, 35)
        assert not sdss_astrombad(77, 6, 77)
        assert sdss_astrombad(85, 1, 15)
        assert (sdss_astrombad(np.array([77, 85, 251]),
                np.array([1, 2, 3]), np.array([20, 15, 151])) ==
                np.array([False, True, True])).all()

    def test_sdss_astrombad_raises(self):
        with raises(ValueError):
            foo = sdss_astrombad(77, 32, 20)
        with raises(ValueError):
            foo = sdss_astrombad(-1, 1, 20)
        with raises(ValueError):
            foo = sdss_astrombad(2**17, 1, 20)
        with raises(ValueError):
            foo = sdss_astrombad(-2, 1, 20)
        with raises(ValueError):
            foo = sdss_astrombad(251, 1, 2**16)
        with raises(ValueError):
            foo = sdss_astrombad(np.array([77, 85, 251]), np.array([1]),
                                np.array([20, 15, 151]))
        with raises(ValueError):
            foo = sdss_astrombad(np.array([77, 85, 251]), np.array([1, 2, 3]),
                                np.array([20]))

    @remote_data
    def test_sdss_astrombad_remote(self):
        pydl.pydlutils.sdss.opbadfields = None
        assert not sdss_astrombad(77, 1, 20)
        assert sdss_astrombad(77, 3, 35)
        assert not sdss_astrombad(77, 6, 77)

    def test_sdss_objid(self):
        assert sdss_objid(3704, 3, 91, 146) == 1237661382772195474
        run = np.array([3704, 1000], dtype=np.int64)
        camcol = np.array([3, 6], dtype=np.int64)
        field = np.array([91, 77], dtype=np.int64)
        obj = np.array([146, 123], dtype=np.int64)
        assert (np.array([1237661382772195474, 1237649770790322299],
                         dtype=np.int64) ==
                sdss_objid(run, camcol, field, obj)).all()
        #
        # Exceptions
        #
        with raises(ValueError):
            objid = sdss_objid(run, 3, 91, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, camcol, 91, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, field, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, obj)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, 146, rerun=np.array([137, 301]))
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, 146, skyversion=np.array([2, 3]))
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, 146, skyversion=-2)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, 146, skyversion=16)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, 146, rerun=-2)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, 146, rerun=2**11)
        with raises(ValueError):
            objid = sdss_objid(-2, 3, 91, 146)
        with raises(ValueError):
            objid = sdss_objid(2**16, 3, 91, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, 0, 91, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, 7, 91, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, -2, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 2**12, 146)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, -2)
        with raises(ValueError):
            objid = sdss_objid(3704, 3, 91, 2**16)

    def test_sdss_specobjid(self):
        s = sdss_specobjid(4055, 408, 55359, 'v5_7_0')
        assert s == 4565636362342690816
        s = sdss_specobjid(4055, 408, 55359, '700')
        assert s == 4565636362342690816
        s = sdss_specobjid(4055, 408, 55359, 700)
        assert s == 4565636362342690816
        s = sdss_specobjid(4055, 408, 55359, 'v5_7_0', line=137)
        assert s == 4565636362342690816 + 137
        s = sdss_specobjid(4055, 408, 55359, 'v5_7_0', index=137)
        assert s == 4565636362342690816 + 137
        #
        # Because 5-digit plate IDs...
        #
        s = sdss_specobjid(10000, 1, 57346, 'v5_10_0')
        assert s == 11258999466550599680
        #
        # Exceptions
        #
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, 'v5_7_0', line=137,
                               index=137)
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, 'V5.7.0')
        with raises(ValueError):
            s = sdss_specobjid(-1, 408, 55359, 'v5_7_0')
        with raises(ValueError):
            s = sdss_specobjid(2**15, 408, 55359, 'v5_7_0')
        with raises(ValueError):
            s = sdss_specobjid(4055, -1, 55359, 'v5_7_0')
        with raises(ValueError):
            s = sdss_specobjid(4055, 2**12, 55359, 'v5_7_0')
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 49999, 'v5_7_0')
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 2**15, 'v5_7_0')
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, -1)
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, 2**15)
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, 'v5_7_0', line=-1)
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, 'v5_7_0', line=2**10)
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, 'v5_7_0', index=-1)
        with raises(ValueError):
            s = sdss_specobjid(4055, 408, 55359, 'v5_7_0', index=2**10)

    def test_unwrap_specobjid(self):
        s = 4565636362342690816
        d = unwrap_specobjid(np.array([s], dtype=np.uint64))
        assert d.plate == 4055
        assert d.mjd == 55359
        assert d.fiber == 408
        assert d.run2d == 'v5_7_0'
        assert d.line == 0
        d = unwrap_specobjid(np.array(['4565636362342690816']))
        assert d.plate == 4055
        assert d.mjd == 55359
        assert d.fiber == 408
        assert d.run2d == 'v5_7_0'
        assert d.line == 0
        si = 4565636362342690816 + 137
        d = unwrap_specobjid(np.array([si], dtype=np.uint64))
        assert d.plate == 4055
        assert d.mjd == 55359
        assert d.fiber == 408
        assert d.run2d == 'v5_7_0'
        assert d.line == 137
        d = unwrap_specobjid(np.array([si], dtype=np.uint64),
                             specLineIndex=True)
        assert d.plate == 4055
        assert d.mjd == 55359
        assert d.fiber == 408
        assert d.run2d == 'v5_7_0'
        assert d.index == 137
        d = unwrap_specobjid(np.array([si], dtype=np.uint64),
                             run2d_integer=True)
        assert d.plate == 4055
        assert d.mjd == 55359
        assert d.fiber == 408
        assert d.run2d == 700
        assert d.line == 137
        #
        # Because 5-digit plate IDs...
        #
        s5 = 11258999466550599680
        d = unwrap_specobjid(np.array([s5], dtype=np.uint64))
        assert d.plate == 10000
        assert d.mjd == 57346
        assert d.fiber == 1
        assert d.run2d == 'v5_10_0'
        assert d.line == 0
        #
        # Exceptions
        #
        with raises(ValueError):
            d = unwrap_specobjid(np.array([4565636362342690816],
                                 dtype=np.int64))
