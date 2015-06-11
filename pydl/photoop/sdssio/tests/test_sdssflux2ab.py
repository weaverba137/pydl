# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_sdssflux2ab():
    import numpy as np
    from .. import sdssflux2ab
    correction = np.array([-0.042, 0.036, 0.015, 0.013, -0.002])
    mags = np.zeros((2,5),dtype='d')
    mags[0,:] = 18.0
    mags[1,:] = 19.0
    ab = sdssflux2ab(mags,magnitude=True)
    assert (ab == (mags + correction)).all()
    flux = 10**((22.5 - mags)/2.5) # nanomaggies
    ab = sdssflux2ab(flux)
    assert np.allclose(ab,10**((22.5 - (mags + correction))/2.5))
    ivar = 1.0/flux
    ab = sdssflux2ab(ivar,ivar=True)
    assert np.allclose(ab,ivar/(10**(-2.0*correction/2.5)))
