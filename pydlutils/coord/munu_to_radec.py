#
#
#
def munu_to_radec(mu, nu, **kwargs):
    """Convert from SDSS great circle coordinates to equatorial coordinates.
    """
    import numpy as np
    from pydlutils.coord import stripe_to_eta
    from pydlutils.goddard.misc import cirrange
    if 'stripe' in kwargs:
        node = 95.0
        incl = stripe_to_incl(kwargs['stripe'])
    elif 'node' in kwargs and 'incl' in kwargs:
        node = kwargs['node']
        incl = kwargs['incl']
    else:
        raise ValueError('Must specify either STRIPE or NODE,INCL!')
    if mu.size != nu.size:
        raise ValueError('Number of elements in MU and NU must agree!')
    sinnu = np.sin(np.deg2rad(nu))
    cosnu = np.cos(np.deg2rad(nu))
    sini = np.sin(np.deg2rad(incl))
    cosi = np.cos(np.deg2rad(incl))
    sinmu = np.sin(np.deg2rad(mu-node))
    cosmu = np.cos(np.deg2rad(mu-node))
    xx = cosmu * cosnu
    yy = sinmu * cosnu * cosi - sinnu * sini
    zz = sinmu * cosnu * sini + sinnu * cosi
    ra = np.rad2deg(np.arctan2(yy,xx)) + node
    dec = np.rad2deg(np.arcsin(zz))
    if 'phi' in kwargs:
        phi = np.rad2deg(np.arctan2(cosmu * sini,
            (-sinmu * sinnu * sini + cosnu * cosi)*cosnu))
        return (ra,dec,phi)
    else:
        return (ra,dec)

