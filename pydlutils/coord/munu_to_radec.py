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
    sinnu = np.sin(np.radians(nu))
    cosnu = np.cos(np.radians(nu))
    sini = np.sin(np.radians(incl))
    cosi = np.cos(np.radians(incl))
    sinmu = np.sin(np.radians(mu-node))
    cosmu = np.cos(np.radians(mu-node))
    xx = cosmu * cosnu
    yy = sinmu * cosnu * cosi - sinnu * sini
    zz = sinmu * cosnu * sini + sinnu * cosi
    ra = np.degrees(np.arctan2(yy,xx)) + node
    dec = np.degrees(np.arcsin(zz))
    ra = np.zeros(mu.shape,dtype=mu.dtype)
    dec = np.zeros(nu.shape,dtype=mu.dtype)
    if 'phi' in kwargs:
        phi = np.degrees(np.arctan2(cosmu * sini,
            (-sinmu * sinnu * sini + cosnu * cosi)*cosnu))
        return (ra,dec,phi)
    else:
        return (ra,dec)

