# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_spheregroup():
    import numpy as np
    from .. import spheregroup
    np.random.seed(137)
    Ngroup = 3
    N = 50
    x0 = np.concatenate( ( np.random.normal(loc=1,scale=0.05,size=(N,)),
                           np.random.normal(loc=-1,scale=0.05,size=(N,)),
                           np.random.normal(loc=1,scale=0.05,size=(N,)),
                         )).reshape((N*Ngroup,))
    y0 = np.concatenate( ( np.random.normal(loc=1,scale=0.05,size=(N,)),
                           np.random.normal(loc=-1,scale=0.05,size=(N,)),
                           np.random.normal(loc=-1,scale=0.05,size=(N,)),
                         )).reshape((N*Ngroup,))
    z0 = np.concatenate( ( np.random.normal(loc=1,scale=0.05,size=(N,)),
                           np.random.normal(loc=-1,scale=0.05,size=(N,)),
                           np.random.normal(loc=0,scale=0.05,size=(N,)),
                         )).reshape((N*Ngroup,))
    foo = np.arange(N*Ngroup)
    np.random.shuffle(foo)
    x = x0[foo]
    y = y0[foo]
    z = z0[foo]
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.degrees(np.arccos(z/r))
    phi = np.degrees(np.arctan2(y,x))
    ra = np.where(phi < 0, phi + 360.0,phi)
    dec = 90.0 - theta
    # group = spheregroup(ra,dec,5.0)
