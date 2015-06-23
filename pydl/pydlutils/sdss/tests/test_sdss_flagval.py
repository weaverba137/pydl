# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_sdss_flagval():
    from os.path import dirname, join
    from astropy.tests.helper import raises
    import pydl.pydlutils.sdss
    from ..set_maskbits import set_maskbits
    pydl.pydlutils.sdss.maskbits = set_maskbits(maskbits_file=join(dirname(__file__),'t','testMaskbits.par'))
    from .. import sdss_flagval
    val = sdss_flagval('TARGET','ROSAT_A')
    assert val == 2**9
    val = sdss_flagval('ANCILLARY_TARGET1',['BLAZGX','ELG','BRIGHTGAL'])
    assert val == 2310346608843161600
    with raises(KeyError):
        val = sdss_flagval('TARGET','ROSAT_Q')
    with raises(KeyError):
        val = sdss_flagval('ABADMASK',"ABADFLAG")
