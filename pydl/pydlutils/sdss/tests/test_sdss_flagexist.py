# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_sdss_flagexist():
    from os.path import dirname, join
    import pydl.pydlutils.sdss
    from ..set_maskbits import set_maskbits
    pydl.pydlutils.sdss.maskbits = set_maskbits(maskbits_file=join(dirname(__file__),'t','testMaskbits.par'))
    from .. import sdss_flagexist
    assert sdss_flagexist('TARGET','ROSAT_A')
    assert sdss_flagexist('TARGET',['rosat_a','rosat_b','rosat_c','rosat_d'])
    l, f = sdss_flagexist('target','FOOBAR',flagexist=True)
    assert not l
    assert f
    l, which = sdss_flagexist('TARGET',['rosat_a','rosat_b','rosat_c','rosat_d','foobar'],whichexist=True)
    assert not l
    assert tuple(which) == (True,True,True,True,False)
    l, f, which = sdss_flagexist('TARGET',['rosat_a','rosat_b','rosat_c','rosat_d','foobar'],
        flagexist=True, whichexist=True)
    assert not l
    assert f
    assert tuple(which) == (True,True,True,True,False)
