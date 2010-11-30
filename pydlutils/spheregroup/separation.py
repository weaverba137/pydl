#
# $Id$
#
def separation(ra1,dec1,ra2,dec2):
    """Takes two points on the celestial sphere IN DECIMAL DEGREES and returns
    their separation IN DECIMAL DEGREES.
    """
    from math import degrees, radians, sin, cos, acos, sqrt
    lng1 = radians(ra1)
    lat1 = radians(dec1)
    lng2 = radians(ra2)
    lat2 = radians(dec2)
    if abs(lat1-lat2) + abs(lng1-lng2) > 1.0e-3:
        sep = acos(cos(lat1)*cos(lat2)*cos(lng1-lng2) + sin(lat1)*sin(lat2))
    else:
        sep = sqrt(cos(lat1)*cos(lat2)*(lng1-lng2)*(lng1-lng2)
            + (lat1-lat2)*(lat1-lat2))
    return degrees(sep)
