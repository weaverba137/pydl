# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from astropy.tests.helper import raises
from ..mangle import is_cap_used

class TestMangle(object):
    """Test the functions in pydl.pydlutils.mangle.
    """

    def setup(self):
        pass

    def teardown(self):
        pass

    def test_is_cap_used(self):
        assert is_cap_used(1 << 2, 2)
        assert not is_cap_used(1 << 2, 1)
