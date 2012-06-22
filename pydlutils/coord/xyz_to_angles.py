#
#
#
def xyz_to_angles(x,y,z):
    """Convert Cartesion coordinates to spherical coordinates (r,phi,theta)

    This uses the standard mathematical definition of spherical coordinates.
    To convert theta to Dec, use Dec = 90 - theta.
    """
    import numpy as np
    r = np.sqrt(x*x + y*y + z*z)
    theta = np.rad2deg(np.arccos(z/r))
    phi = np.rad2deg(np.arctan2(y,x))
    w = phi < 0
    if w.sum() > 0:
        phi[w] += 360.0
    return (r,phi,theta)

