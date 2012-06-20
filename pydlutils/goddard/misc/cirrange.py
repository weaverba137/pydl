#
#
#
def cirrange(ang,radians=False):
    from math import pi
    if radians:
        cnst = pi * 2.0
    else:
        cnst = 360.0
    #
    # The modulo operator automatically deals with negative values
    #
    return ang % cnst

