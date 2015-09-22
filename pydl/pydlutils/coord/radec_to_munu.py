# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from astropy.coordinates import Angle, frame_transform_graph, FunctionTransform, ICRS
from . import SDSSMuNu
#
@frame_transform_graph.transform(FunctionTransform, ICRS, SDSSMuNu)
def radec_to_munu(icrs_frame,munu):
    """Convert from equatorial coordinates to SDSS great circle coordinates.

    Parameters
    ----------
    icrs_frame : :class:`~astropy.coordinates.ICRS`
        Equatorial coordinates (RA, Dec).

    Returns
    -------
    radec_to_munu : :class:`~pydl.pydlutils.coord.SDSSMuNu`
        SDSS great circle coordinates (mu, nu).
    """
    from astropy import units as u
    import numpy as np
    # from pydlutils.coord import stripe_to_eta
    # from pydlutils.goddard.misc import cirrange
    # if 'stripe' in kwargs:
    #     node = 95.0
    #     incl = stripe_to_incl(kwargs['stripe'])
    # elif 'node' in kwargs and 'incl' in kwargs:
    #     node = kwargs['node']
    #     incl = kwargs['incl']
    # else:
    #     raise ValueError('Must specify either STRIPE or NODE,INCL!')
    # if ra.size != dec.size:
    #     raise ValueError('Number of elements in RA and DEC must agree!')
    sinra = np.sin((icrs_frame.ra - munu.node).to(u.radian).value)
    cosra = np.cos((icrs_frame.ra - munu.node).to(u.radian).value)
    sindec = np.sin(icrs_frame.dec.to(u.radian).value)
    cosdec = np.cos(icrs_frame.dec.to(u.radian).value)
    sini = np.sin(munu.incl.to(u.radian).value)
    cosi = np.cos(munu.incl.to(u.radian).value)
    x1 = cosdec * cosra
    y1 = cosdec * sinra
    z1 = sindec
    x2 = x1
    y2 = y1 * cosi + z1 * sini
    z2 = -y1 * sini + z1 * cosi
    mu = Angle(np.arctan2(y2,x2),unit=u.radian) + munu.node
    nu = Angle(np.arcsin(z2),unit=u.radian)
    # if 'phi' in kwargs:
    #     sinnu = np.sin(np.deg2rad(nu))
    #     cosnu = np.cos(np.deg2rad(nu))
    #     sinmu = np.sin(np.deg2rad(mu-node))
    #     cosmu = np.cos(np.deg2rad(mu-node))
    #     phi = np.rad2deg(np.arctan2(cosmu * sini,
    #         (-sinmu * sinnu * sini + cosnu * cosi)*cosnu))
    #     return (ra,dec,phi)
    # else:
    #     return (ra,dec)
    return SDSSMuNu(mu=mu,nu=nu,stripe=munu.stripe)
