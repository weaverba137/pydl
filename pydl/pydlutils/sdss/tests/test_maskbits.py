# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.tests.helper import raises
from os.path import dirname, join
import pydl.pydlutils.sdss
from ..set_maskbits import set_maskbits
from .. import sdss_flagexist, sdss_flagname, sdss_flagval


class TestMaskbits(object):
    """Test all the functions related to bitmasks.
    """

    def setup(self):
        self.data_dir = join(dirname(__file__), 't')
        pydl.pydlutils.sdss.maskbits = set_maskbits(
            maskbits_file=join(self.data_dir, 'testMaskbits.par'))

    def teardown(self):
        pass

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
