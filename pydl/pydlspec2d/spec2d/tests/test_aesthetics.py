# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
import pytest
#
def test_aesthetics():
    from .. import aesthetics
    from ... import Pydlspec2dException
    import numpy as np
    np.random.seed(137)
    flux = np.random.rand(100)
    ivar = np.random.rand(100)
    #
    # No bad
    #
    f = aesthetics(flux,ivar)
    assert (f == flux).all()
    #
    # Bad points
    #
    ivar[ivar < 0.1] = 0.0
    #
    # Bad method
    #
    with pytest.raises(Pydlspec2dException):
        f = aesthetics(flux,ivar,'badmethod')
    #
    # Nothing
    #
    f = aesthetics(flux,ivar,'nothing')
    assert (f == flux).all()
