#
#
#
def angles_to_xyz(r,phi,theta):
    """Convert spherical coordinates (r,phi,theta) into Cartesion coordinates

    This uses the standard mathematical definition of spherical coordinates.
    To convert e.g. RA, Dec or longitude, latitude, pass RA, 90 - Dec.
    """
    import numpy as np
    x = r * np.cos(np.deg2rad(phi)) * np.sin(np.deg2rad(theta))
    y = r * np.sin(np.deg2rad(phi)) * np.sin(np.deg2rad(theta))
    z = r * np.cos(np.deg2rad(theta))
    return (x,y,z)

