# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_default_skyversion():
    from .. import default_skyversion
    assert default_skyversion() == 2
