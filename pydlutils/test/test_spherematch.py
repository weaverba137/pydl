#
# $Id$
#
def test_spherematch(save=False):
    import numpy as np
    import os
    import os.path
    import pyfits
    import time
    from pydlutils.spheregroup import spherematch
    from pydlutils.misc import djs_readcol
    searchrad = 3.0/3600.0
    n = 20
    ra1 = 360.0*np.random.random((n,))
    dec1 = 90.0 - np.rad2deg(np.arccos(2.0*np.random.random((n,)) - 1.0))
    ra2 = ra1 + np.random.normal(0,1.0/3600.0)
    dec2 = dec1 + np.random.normal(0,1.0/3600.0)
    foo = np.arange(n)
    np.random.shuffle(foo)
    i1, i2, d12 = spherematch(ra1, dec1, ra2[foo], dec2[foo], searchrad, maxmatch=0)
    print ra1
    print dec1
    print ra2[foo]
    print dec2[foo]
    print foo
    print i1
    print i2
    print d12
    if save:
        np.savetxt('file.out',np.vstack((ra1,dec1,ra2[foo],dec2[foo])).T)
    ra, dec = djs_readcol('file.in',format='(D,D)')
    run2d = '26'
    run1d = None
    plates = pyfits.open(os.path.join(os.getenv('SPECTRO_REDUX'),'plates-dr8.fits'))
    plist = plates[1].data
    plates.close()
    qdone = plist.field('STATUS1D') == 'Done'
    qdone2d = plist.field('RUN2D').strip() == run2d
    if run1d is None:
        qdone1d = np.ones(plist.size,dtype='bool')
    else:
        qdone1d = plist.field('RUN1D').strip() == run1d
    qfinal = qdone & qdone2d & qdone1d
    idone = np.arange(plist.size)[qfinal]
    t0 = time.time()
    imatch1, itmp, dist12 = spherematch(ra, dec,
        plist[idone].field('RACEN'), plist[idone].field('DECCEN'),
        searchrad+1.55,maxmatch=0,debug=True)
    t1 = time.time()
    print "Elapsed time = %f s." % (t1-t0)
    foo = plist[idone]
    print imatch1
    print itmp
    print dist12
    t2 = time.time()
    imatch1, itmp, dist12 = spherematch(
        plist[idone].field('RACEN'), plist[idone].field('DECCEN'), ra, dec,
        searchrad+1.55,maxmatch=0)
    t3 = time.time()
    print "Elapsed time = %f s." % (t3-t2)
    print imatch1
    print itmp
    print dist12

    return
