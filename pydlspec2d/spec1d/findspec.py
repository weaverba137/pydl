#
# $Id$
#
def findspec(*args,**kwargs):
    """Find SDSS/BOSS spectra that match a given RA, Dec
    """
    import os
    import os.path
    import glob
    import pyfits
    import numpy as np
    from pydlutils.misc import djs_readcol, struct_print
    from pydlspec2d import Pydlspec2dException
    import pydlspec2d.spec1d # get the findspec_cache dictionary
    #
    # Set up default values
    #
    if 'sdss' in kwargs:
        topdir = os.getenv('SPECTRO_REDUX')
        run2d = '26'
        run1d = None
    else:
        if 'topdir' in kwargs:
            topdir = kwargs['topdir']
        else:
            topdir = os.getenv('BOSS_SPECTRO_REDUX')
        if 'run2d' in kwargs:
            run2d = str(kwargs['run2d'])
        else:
            run2d = os.getenv('RUN2D')
        if 'run1d' in kwargs:
            run1d = str(kwargs['run1d'])
        else:
            run1d = os.getenv('RUN1D')
    if pydlspec2d.spec1d.findspec_cache is None:
        pydlspec2d.spec1d.findspec_cache = { 'lasttopdir':topdir,
            'plist':None }
    if (pydlspec2d.spec1d.findspec_cache['plist'] is None or
        topdir != pydlspec2d.spec1d.findspec_cache['lasttopdir']):
        pydlspec2d.spec1d.findspec_cache['lasttopdir'] = topdir
        platelist_file = "%s/platelist.fits" % topdir
        plates_files = glob.glob("%s/plates-*.fits" % topdir)
        plist = None
        if os.path.exists(platelist_file):
            platelist = pyfits.open(platelist_file)
            plist = platelist[1].data
            platelist.close()
        if os.path.exists(plates_files[0]):
            plates = pyfits.open(plates_files[0])
            plist = plates[1].data
            plates.close()
        if plist is None:
            raise Pydlspec2dException("Plate list (platelist.fits or plates-*.fits) not found in %s." % topdir)
        else:
            pydlspec2d.spec1d.findspec_cache['plist'] = plist
    qdone = plist.field('STATUS1D') == 'Done'
    qdone2d = plist.field('RUN2D').strip() == run2d
    if run1d is None:
        qdone1d = np.ones(plist.size,dtype='bool')
    else:
        qdone1d = plist.field('RUN1D').strip() == run1d
    qfinal = qdone & qdone2d & qdone1d
    if not qfinal.any():
        print "No reduced plates!"
        return None
    #
    # If there are positional arguments, interpret these as RA, Dec
    #
    if len(args) == 2:
        ra = args[0]
        dec = args[1]
    #
    # Read RA, Dec from infile if set
    #
    if 'infile' in kwargs:
        ra, dec = djs_readcol(kwargs['infile'],format='(D,D)')
    if 'searchrad' in kwargs:
        searchrad = float(kwargs['searchrad'])
    else:
        searchrad = 3.0/3600.0
    #
    # Create output structure
    #
    slist = np.zeros(plist.size,dtype=[('PLATE','i4'),
        ('MJD','i4'),('FIBERID','i4'),
        ('RA','f8'),('DEC','f8'),('MATCHRAD','f4')])
    imatch1, itmp, dist12 = spherematch(ra, dec, plist[qfinal].field('RACEN'), plist[qfinal].field('DECCEN'),searchrad+1.55)
    if (imatch[0] == -1):
        return None
    imatch2 = qfinal[itmp]
    try:
        n_total = plist.field('N_TOTAL')
    except KeyError:
        n_total = np.zeros(plist.size,dtype='i4') + 640
    #
    # Print to terminal or output file
    #
    if 'print' in kwargs:
        foo = struct_print(slist)
    if 'outfile' in kwargs:
        foo = struct_print(slist,filename=outfile)
    return slist
