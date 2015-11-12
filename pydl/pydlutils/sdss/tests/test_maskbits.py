# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.tests.helper import raises
from os.path import dirname, join
import pydl.pydlutils.sdss
from ..set_maskbits import set_maskbits
from .. import sdss_flagname, sdss_flagval


class TestMaskbits(object):
    """Test all the functions related to bitmasks.
    """

    def setup(self):
        self.data_dir = join(dirname(__file__), 't')
        pydl.pydlutils.sdss.maskbits = set_maskbits(
            maskbits_file=join(self.data_dir, 'testMaskbits.par'))

    def teardown(self):
        pass

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
