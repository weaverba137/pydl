#
#
#
def gcirc(ra1,dec1,ra2,dec2,units=2):
    """Computes rigorous great circle arc distances.

    units = 0: everything is already in radians
    units = 1: RA in hours, dec in degrees, distance in arcsec.
    units = 2: RA, dec in degrees, distance in arcsec (default)

    The formula below is the one best suited to handling small angular
    separations.  See:
    http://en.wikipedia.org/wiki/Great-circle_distance
    """
    import numpy as np
    if units == 0:
        rarad1 = ra1
        dcrad1 = dec1
        rarad2 = ra2
        dcrad2 = dec2
    elif units == 1:
        rarad1 = np.deg2rad(15.0*ra1)
        dcrad1 = np.deg2rad(dec1)
        rarad2 = np.deg2rad(15.0*ra2)
        dcrad2 = np.deg2rad(dec2)
    elif units == 2:
        rarad1 = np.deg2rad(ra1)
        dcrad1 = np.deg2rad(dec1)
        rarad2 = np.deg2rad(ra2)
        dcrad2 = np.deg2rad(dec2)
    else:
        raise ValueError('units must be 0, 1 or 2!')
    deldec2 = (dcrad2-dcrad1)/2.0
    delra2 =  (rarad2-rarad1)/2.0
    sindis = np.sqrt( np.sin(deldec2)*np.sin(deldec2) +
        np.cos(dcrad1)*np.cos(dcrad2)*np.sin(delra2)*np.sin(delra2) )
    dis = 2.0*np.arcsin(sindis)
    if units == 0:
        return dis
    else:
        return np.rad2deg(dis)*3600.0

