# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_func_fit():
    from .. import func_fit
    import numpy as np
    from astropy.tests.helper import raises
    np.random.seed(137)
    x = np.linspace(-5,5,50)
    y = x**2 + 2*x + 1 + 0.05*np.random.randn(50)
    #
    # Bad inputs
    #
    with raises(ValueError):
        foo = func_fit(x,np.array([1,2,3]),3)
    with raises(ValueError):
        foo = func_fit(x,y,3,invvar=np.array([0,0,0]))
    with raises(ValueError):
        foo = func_fit(x,y,3,inputfunc=np.array([1.0,1.0,1.0]))
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
    invvar[invvar < 2] = 0
    res,yfit = func_fit(x,y,3,invvar=invvar,function_name='poly')
    # assert np.allclose(res,np.array([ 0.99665423,  1.9945388 ,  1.00172303]))
    assert np.allclose(res,np.array([ 0.99996197,  1.99340315,  1.00148004]))
    #
    # Fixed points
    #
    res,yfit = func_fit(x,y,3,invvar=invvar,function_name='poly',ia=np.array([False,True,True]),inputans=np.array([1.0,0,0]))
    # assert np.allclose(res,np.array([ 1.        ,  1.99454782,  1.00149949]))
    assert np.allclose(res,np.array([ 1.        ,  1.99340359,  1.00147743]))
    res,yfit = func_fit(x,y,3,invvar=invvar,function_name='poly',ia=np.array([False,True,False]),inputans=np.array([1.0,0,1.0]))
    # assert np.allclose(res,np.array([ 1.        ,  1.99403239,  1.        ]))
    assert np.allclose(res,np.array([ 1.        ,  1.99735654,  1.        ]))
    #
    # inputfunc
    #
    res,yfit = func_fit(x,y,3,invvar=invvar,function_name='poly',inputfunc=np.ones(x.shape,dtype=x.dtype))
    # assert np.allclose(res,np.array([ 0.99665423,  1.9945388 ,  1.00172303]))
    assert np.allclose(res,np.array([ 0.99996197,  1.99340315,  1.00148004]))
    #
    # Generate inputans
    #
    y = x**2 + 2*x + 0.05*np.random.randn(50)
    res,yfit = func_fit(x,y,3,invvar=invvar,function_name='poly',ia=np.array([False,True,True]))
    assert np.allclose(res,np.array([ 0.        ,  1.99994188,  0.99915111]))
