# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from __future__ import print_function
#
def findspec(*args,**kwargs):
    """Find SDSS/BOSS spectra that match a given RA, Dec.

    Parameters
    ----------

    Returns
    -------
    """
    import os
    import os.path
    import glob
    from astropy.io import ascii, fits
    import numpy as np
    from ... import uniq
    from ...pydlutils.misc import struct_print
    from ...pydlutils.spheregroup import spherematch
    from .. import Pydlspec2dException
    import pydl.pydlspec2d.spec1d # get the findspec_cache dictionary
    #
    # Set up default values
    #
    if 'sdss' in kwargs:
        topdir = os.getenv('SPECTRO_REDUX')
        run2d = '26'
        run1d = ''
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
    if pydl.pydlspec2d.spec1d.findspec_cache is None:
        pydl.pydlspec2d.spec1d.findspec_cache = { 'lasttopdir':topdir,
            'plist':None }
    if (pydl.pydlspec2d.spec1d.findspec_cache['plist'] is None or
        topdir != pydl.pydlspec2d.spec1d.findspec_cache['lasttopdir']):
        pydl.pydlspec2d.spec1d.findspec_cache['lasttopdir'] = topdir
        platelist_file = os.path.join(topdir,"platelist.fits")
        plates_files = glob.glob(os.path.join(topdir,"plates-*.fits"))
        plist = None
        if os.path.exists(platelist_file):
            platelist = fits.open(platelist_file)
            plist = platelist[1].data
            platelist.close()
        if len(plates_files) > 0:
            plates = fits.open(plates_files[0])
            plist = plates[1].data
            plates.close()
        if plist is None:
            raise Pydlspec2dException("Plate list (platelist.fits or plates-*.fits) not found in {0}.".format(topdir))
        else:
            pydl.pydlspec2d.spec1d.findspec_cache['plist'] = plist
    qdone = plist.field('STATUS1D') == 'Done'
    qdone2d = plist.field('RUN2D').strip() == run2d
    if run1d == '':
        qdone1d = np.ones(plist.size,dtype='bool')
    else:
        qdone1d = plist.field('RUN1D').strip() == run1d
    qfinal = qdone & qdone2d & qdone1d
    if not qfinal.any():
        print("No reduced plates!")
        return None
    idone = np.arange(plist.size)[qfinal]
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
        infile_data = ascii.read(kwargs['infile'],names=['ra','dec'])
        ra = infile_data["ra"].data
        dec = infile_data["dec"].data
    if 'searchrad' in kwargs:
        searchrad = float(kwargs['searchrad'])
    else:
        searchrad = 3.0/3600.0
    #
    # Create output structure
    #
    slist_type = np.dtype([('PLATE','i4'),('MJD','i4'),('FIBERID','i4'),
        ('RA','f8'),('DEC','f8'),('MATCHRAD','f8')])
    #
    # Match all plates with objects
    #
    imatch1, itmp, dist12 = spherematch(ra, dec,
        plist[qfinal].field('RACEN'), plist[qfinal].field('DECCEN'),
        searchrad+1.55,maxmatch=0)
    if imatch1.size == 0:
        return None
    imatch2 = idone[itmp]
    #
    # Read all relevant plates
    #
    try:
        n_total = plist.field('N_TOTAL')
    except KeyError:
        n_total = np.zeros(plist.size,dtype='i4') + 640
    iplate = imatch2[uniq(imatch2,imatch2.argsort())]
    i0 = 0
    plugmap = np.zeros(n_total[iplate].sum(),
        dtype=[('PLATE','i4'),('MJD','i4'),('FIBERID','i4'),
        ('RA','d'),('DEC','d')])
    for i in range(iplate.size):
        spplate = pydl.pydlspec2d.spec1d.readspec(plist[iplate[i]].field('PLATE'),mjd=plist[iplate[i]].field('MJD'),
            topdir=topdir,run2d=run2d,run1d=run1d)
        index_to = i0 + np.arange(n_total[iplate[i]],dtype='i4')
        plugmap['PLATE'][index_to] = plist[iplate[i]].field('PLATE')
        plugmap['MJD'][index_to] = plist[iplate[i]].field('MJD')
        plugmap['FIBERID'][index_to] = spplate['plugmap']['FIBERID']
        plugmap['RA'][index_to] = spplate['plugmap']['RA']
        plugmap['DEC'][index_to] = spplate['plugmap']['DEC']
        i0 += n_total[iplate[i]]
    i1, i2, d12 = spherematch(ra, dec, plugmap['RA'], plugmap['DEC'], searchrad, maxmatch=0)
    if i1.size == 0:
        return None
    if 'best' in kwargs:
        #
        # Return only best match per object
        #
        slist = np.zeros(ra.size,dtype=slist_type)
        spplate = pydl.pydlspec2d.spec1d.readspec(plugmap[i2]['PLATE'],plugmap[i2]['FIBERID'],
            mjd=plugmap[i2]['MJD'],topdir=topdir, run2d=run2d, run1d=run1d)
        sn = spplate['zans']['SN_MEDIAN']
        isort = (i1 + np.where(sn > 0, sn, 0)/(sn+1.0).max()).argsort()
        i1 = i1[isort]
        i2 = i2[isort]
        d12 = d12[isort]
        iuniq = uniq(i1)
        slist[i1[iuniq]]['PLATE'] = plugmap[i2[iuniq]]['PLATE']
        slist[i1[iuniq]]['MJD'] = plugmap[i2[iuniq]]['MJD']
        slist[i1[iuniq]]['FIBERID'] = plugmap[i2[iuniq]]['FIBERID']
        slist[i1[iuniq]]['RA'] = plugmap[i2[iuniq]]['RA']
        slist[i1[iuniq]]['DEC'] = plugmap[i2[iuniq]]['DEC']
        slist[i1[iuniq]]['MATCHRAD'] = d12[iuniq]
    else:
        #
        # Return all matches
        #
        slist = np.zeros(i1.size,dtype=slist_type)
        slist['PLATE'] = plugmap[i2]['PLATE']
        slist['MJD'] = plugmap[i2]['MJD']
        slist['FIBERID'] = plugmap[i2]['FIBERID']
        slist['RA'] = plugmap[i2]['RA']
        slist['DEC'] = plugmap[i2]['DEC']
        slist['MATCHRAD'] = d12
    #
    # Print to terminal or output file
    #
    if 'print' in kwargs:
        foo = struct_print(slist)
    if 'outfile' in kwargs:
        foo = struct_print(slist,filename=outfile)
    return slist
