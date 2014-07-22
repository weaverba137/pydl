# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def test_struct_print():
    import numpy as np
    from .. import struct_print
    n = 20
    slist = np.zeros(n,dtype=[('PLATE','i4'),('MJD','i4'),('FIBERID','i4'),('RA','f8'),('DEC','f8'),('MATCHRAD','f4'),('RERUN','S8')])
    slist['RERUN'][0:int(n/2)] = 'v1_2_3'
    slist['RERUN'][int(n/2):n] = 'v1_3_34'
    slist['PLATE'] = np.random.random_integers(10,10000,(n,))
    slist['MJD'] = np.random.random_integers(51000,56000,(n,))
    slist['FIBERID'] = np.random.random_integers(1,1000,(n,))
    slist['RA'] = 360.0*np.random.random((n,))
    slist['DEC'] = 90.0 - np.rad2deg(np.arccos(2.0*np.random.random((n,)) - 1.0))
    slist['MATCHRAD'] = np.random.random((n,))
    #print(slist)
    #lines,css = struct_print(slist,debug=True)
    #print(lines)
    #print(css)
    #lines,css = struct_print(slist,debug=True,html=True)
    #print(lines)
    #print(css)
    assert slist.size == n
