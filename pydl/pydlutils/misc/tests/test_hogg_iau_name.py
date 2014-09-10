# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_hogg_iau_name():
    from numpy import array
    from .. import hogg_iau_name
    assert hogg_iau_name(354.120375,-0.544777778) == 'SDSS J233628.89-003241.2'
    assert hogg_iau_name(354.120375,-0.544777778,prefix='2MASS') == '2MASS J233628.89-003241.2'
    assert hogg_iau_name(354.120375,-0.544777778,prefix='') == 'J233628.89-003241.2'
    assert hogg_iau_name(354.120375,-0.544777778,precision=0) == 'SDSS J233628.8-003241'
    assert hogg_iau_name(354.120375,-0.544777778,precision=2) == 'SDSS J233628.890-003241.20'
    ra = array([354.120375,7.89439,36.31915,110.44730])
    dec = array([-0.544777778,-0.35157,0.47505,39.35352])
    names = hogg_iau_name(ra,dec)
    assert tuple(names) == ('SDSS J233628.89-003241.2',
        'SDSS J003134.65-002105.6','SDSS J022516.59+002830.1','SDSS J072147.35+392112.6')
