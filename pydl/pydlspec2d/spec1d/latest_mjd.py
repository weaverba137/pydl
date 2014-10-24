# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.extern import six

def latest_mjd(plate,**kwargs):
    """Find the most recent MJD associated with a plate.

    Parameters
    ----------
    plate : int or ndarray
        The plate(s) to examine.

    Returns
    -------
    latest_mjd : ndarray
        An array of MJD values for each plate
    """
    import glob
    import re
    import numpy as np
    from . import spec_path
    if isinstance(plate, six.integer_types) or plate.shape == ():
        platevec = np.array([plate],dtype='i4')
    else:
        platevec = plate
    mjd = np.zeros(len(platevec),dtype='i4')
    mjdre = re.compile(r'spPlate-[0-9]{4}-([0-9]{5}).fits')
    unique_plates = np.unique(platevec)
    paths = spec_path(unique_plates,**kwargs)
    for p,q in zip(paths,unique_plates):
        plateglob = "{0}/spPlate-{1:04d}-*.fits".format(p,q)
        bigmjd = 0
        for f in glob.glob(plateglob):
            thismjd = int(mjdre.search(f).groups()[0])
            if thismjd > bigmjd:
                bigmjd = thismjd
        mjd[platevec == q] = bigmjd
    return mjd
