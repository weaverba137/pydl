# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_hogg_iau_name():
    from .. import hogg_iau_name
    assert hogg_iau_name(354.120375,-0.544777778) == 'SDSS J233628.89-003241.2'
