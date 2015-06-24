# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_func_fit():
    from .. import func_fit
    import numpy as np
    from astropy.tests.helper import raises
    np.random.seed(424242)
    x = np.linspace(-5,5,50,dtype='d')
    y = x**2 + 2*x + 1 + 0.05*np.random.randn(50)
    #
    # Bad inputs
    #
    with raises(ValueError):
        foo = func_fit(x,np.array([1,2,3]),3)
    with raises(ValueError):
        foo = func_fit(x,y,3,invvar=np.array([0,0,0]))
    with raises(KeyError):
        foo = func_fit(x,y,3,function_name='npoly')
    #
    # No good points
    #
    invvar = np.zeros(x.shape,dtype=x.dtype)
    res,yfit = func_fit(x,y,3,invvar=invvar)
    assert (res == np.zeros((3,),dtype=x.dtype)).all()
    assert (yfit == 0*x).all()
    #
    # One good point
    #
    invvar[2] = 1.0
    res,yfit = func_fit(x,y,3,invvar=invvar)
    assert (invvar > 0).nonzero()[0][0] == 2
    assert res[0] == y[2]
    assert (yfit == y[2]).all()
    #
    # Various points
    #
    invvar = 1.0/(np.random.random(x.shape)**2)
    assert (invvar < 2).any()
    invvar[invvar < 2] == 0
    res,yfit = func_fit(x,y,3,invvar=invvar,function_name='poly')
    print(res,yfit)
