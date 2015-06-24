# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_func_fit():
    from .. import func_fit
    import numpy as np
    np.random.seed(424242)
    x = np.linspace(-5,5,50,dtype='d')
    y = x**2 + 2*x + 1
    # invvar = 1.0/(np.random.random(x.shape)**2)
    # invvar[invar < 2] == 0
    #
    # No good points
    #
    invvar = np.zeros(x.shape,dtype=x.dtype)
    res,yfit = func_fit(x,y,3,invvar=invvar)
    assert (res == np.zeros((3,),dtype=x.dtype)).all()
    assert (yfit == 0*x).all()
