#
#
#
def latest_mjd(plate,**kwargs):
    """Find the most recent MJD associated with a plate
    """
    import glob
    import re
    import numpy as np
    from pydlspec2d.spec1d import spec_path
    if isinstance(plate,int):
        platevec = np.zeros(1,dtype='i4') + plate
    else:
        platevec = plate
    mjd = np.zeros(len(platevec),dtype='i4')
    mjdre = re.compile(r'spPlate-[0-9]{4}-([0-9]{5}).fits')
    unique_plates = np.unique(platevec)
    paths = spec_path(unique_plates,**kwargs)
    for p,q in zip(paths,unique_plates):
        plateglob = "%s/spPlate-%04d-*.fits" % (p,q)
        bigmjd = 0
        for f in glob.glob(plateglob):
            thismjd = int(mjdre.search(f).groups()[0])
            if thismjd > bigmjd:
                bigmjd = thismjd
        mjd[platevec == q] = bigmjd
    return mjd

