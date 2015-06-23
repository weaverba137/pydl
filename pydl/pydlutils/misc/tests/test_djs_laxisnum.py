# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_djs_laxisnum():
    import numpy as np
    from astropy.tests.helper import raises
    from .. import djs_laxisnum
    #
    # 1d
    #
    assert (np.zeros((4,),dtype='i4') == djs_laxisnum((4,))).all()
    #
    # 2d
    #
    l = np.array([[0,0,0,0],[1,1,1,1],[2,2,2,2],[3,3,3,3]],dtype='i4')
    assert (l == djs_laxisnum((4,4))).all()
    assert (l.T == djs_laxisnum((4,4),iaxis=1)).all()
    with raises(ValueError):
        foo = djs_laxisnum((4,4),iaxis=2)
    #
    # 3d
    #
    l = np.zeros((4,4,4),dtype='i4')
    l[1,:,:] = 1
    l[2,:,:] = 2
    l[3,:,:] = 3
    assert (l == djs_laxisnum((4,4,4))).all()
    assert (l.swapaxes(0,1) == djs_laxisnum((4,4,4),iaxis=1)).all()
    assert (l.swapaxes(0,2) == djs_laxisnum((4,4,4),iaxis=2)).all()
    with raises(ValueError):
        foo = djs_laxisnum((4,4,4),iaxis=3)
    #
    # More d
    #
    with raises(ValueError):
        foo = djs_laxisnum((4,4,4,4))
