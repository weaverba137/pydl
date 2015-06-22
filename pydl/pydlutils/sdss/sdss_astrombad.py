# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdss_astrombad(run,camcol,field,photolog_version='dr10'):
    """For a list of RUN, CAMCOL, FIELD, return whether each field has bad astrometry.

    Parameters
    ----------
    run, camcol, field : int or array of int
        Run, camcol and field.  If arrays are passed,
        all must have the same length.
    photolog_version : str, optional
        Use this version of photolog to obtain the obBadfields.par file,
        if $PHOTOLOG_DIR is not set.

    Returns
    -------
    sdss_astrombad : ndarray of bool
        Array of bool.  ``True`` indicates the field is bad.

    Raises
    ------
    ValueError
        If the sizes of the arrays don't match or if the array values are out of bounds.

    Notes
    -----
    Reads data from ``$PHOTOLOG_DIR/opfiles/opBadFields.par``.

    If there is a problem with one camcol, we assume a
    problem with all camcols.

    """
    from . import opbadfields
    from ..yanny import yanny
    from numpy import array, bool, int64, zeros
    from astropy.utils.data import download_file
    from os import getenv
    from os.path import join
    #
    # Check inputs
    #
    if isinstance(run,int):
        #
        # Assume all inputs are integers & promote to arrays.
        #
        run = array([run],dtype=int64)
        camcol = array([camcol],dtype=int64)
        field = array([field],dtype=int64)
    else:
        #
        # Check that all inputs have the same shape.
        #
        if run.shape != camcol.shape:
            raise ValueError("camcol.shape does not match run.shape!")
        if run.shape != field.shape:
            raise ValueError("field.shape does not match run.shape!")
    #
    # Check ranges of parameters
    #
    if ((run < 0) | (run >= 2**16)).any():
        raise ValueError("run values are out-of-bounds!")
    if ((camcol < 1) | (camcol > 6)).any():
        raise ValueError("camcol values are out-of-bounds!")
    if ((field < 0) | (field >= 2**12)).any():
        raise ValueError("camcol values are out-of-bounds!")
    #
    # Read the file
    #
    if opbadfields is None: # pragma: no cover
        if getenv('PHOTOLOG_DIR') is None:
            if photolog_version == 'trunk' or photolog_version.startswith('branches/'):
                iversion = photolog_version
            else:
                iversion = 'tags/'+photolog_version
            baseurl = 'http://www.sdss3.org/svn/repo/photolog/{0}/opfiles/opBadfields.par'.format(iversion)
            filename = download_file(baseurl,cache=True)
        else:
            filename = join(getenv('PHOTOLOG_DIR'),'opfiles','opBadfields.par')
        astrombadfile = yanny(filename,np=True)
        w = ((astrombadfile['BADFIELDS']['problem'] == 'astrom'.encode()) |
            (astrombadfile['BADFIELDS']['problem'] == 'rotator'.encode()))
        opbadfields = astrombadfile['BADFIELDS'][w]
    #
    # opbadfields already has astrom problems selected at this point
    #
    bad = zeros(run.shape,dtype=bool)
    for row in opbadfields:
        w = ((run == row['run']) &
        (field >= row['firstfield']) & (field < row['lastfield']))
        if w.any():
            bad[w] = True
    return bad
