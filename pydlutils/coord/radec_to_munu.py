#
#
#
def radec_to_munu(ra, dec, **kwargs):
    """Convert from equatorial coordinates to SDSS great circle coordinates.
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
    if ra.size != dec.size:
        raise ValueError('Number of elements in RA and DEC must agree!')
    sinra = np.sin(np.deg2rad(ra-node))
    cosra = np.cos(np.deg2rad(ra-node))
    sindec = np.sin(np.deg2rad(dec))
    cosdec = np.cos(np.deg2rad(dec))
    sini = np.sin(np.deg2rad(incl))
    cosi = np.cos(np.deg2rad(incl))
    x1 = cosdec * cosra
    y1 = cosdec * sinra
    z1 = sindec
    x2 = x1
    y2 = y1 * cosi + z1 * sini
    z2 = -y1 * sini + z1 * cosi
    mu = cirrange(r2d * np.arctan2(y2,x2) + node)
    nu = r2d * np.arcsin(z2)
    if 'phi' in kwargs:
        sinnu = np.sin(np.deg2rad(nu))
        cosnu = np.cos(np.deg2rad(nu))
        sinmu = np.sin(np.deg2rad(mu-node))
        cosmu = np.cos(np.deg2rad(mu-node))
        phi = np.rad2deg(np.arctan2(cosmu * sini,
            (-sinmu * sinnu * sini + cosnu * cosi)*cosnu))
        return (ra,dec,phi)
    else:
        return (ra,dec)

