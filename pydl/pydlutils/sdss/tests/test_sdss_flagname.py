# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_sdss_flagexist():
    from os.path import dirname, join
    from astropy.tests.helper import raises
    import pydl.pydlutils.sdss
    from ..set_maskbits import set_maskbits
    pydl.pydlutils.sdss.maskbits = set_maskbits(maskbits_file=join(dirname(__file__),'t','testMaskbits.par'))
    from .. import sdss_flagname
    names = sdss_flagname('ANCILLARY_TARGET1',2310346608843161600)
    assert tuple(names) == ('BRIGHTGAL', 'BLAZGX', 'ELG')
    names = sdss_flagname('ANCILLARY_TARGET1',2310346608843161600,concat=True)
    assert names == 'BRIGHTGAL BLAZGX ELG'
    with raises(KeyError):
        names = sdss_flagname('ABADMASK',123456789)
