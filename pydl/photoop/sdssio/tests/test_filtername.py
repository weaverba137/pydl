# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_filtername():
    from .. import filtername
    assert filtername(0) == 'u'
    assert filtername(1) == 'g'
    assert filtername(2) == 'r'
    assert filtername(3) == 'i'
    assert filtername(4) == 'z'
    #
    # filtername should return its argument if it's not
    # integer-like
    #
    assert filtername('r') == 'r'
