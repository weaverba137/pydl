# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from astropy.coordinates import Angle, frame_transform_graph, FunctionTransform, ICRS
from . import SDSSMuNu
#
@frame_transform_graph.transform(FunctionTransform, SDSSMuNu, ICRS)
def munu_to_radec(munu,icrs_frame):
    """Convert from SDSS great circle coordinates to equatorial coordinates.

    Parameters
    ----------
    munu : :class:`~pydl.pydlutils.coord.SDSSMuNu`
        SDSS great circle coordinates (mu, nu).

    Returns
    -------
    munu_to_radec : :class:`~astropy.coordinates.ICRS`
        Equatorial coordinates (RA, Dec).
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
    # if mu.size != nu.size:
    #     raise ValueError('Number of elements in MU and NU must agree!')
    sinnu = np.sin(munu.nu.to(u.radian).value)
    cosnu = np.cos(munu.nu.to(u.radian).value)
    sini = np.sin(munu.incl.to(u.radian).value)
    cosi = np.cos(munu.incl.to(u.radian).value)
    sinmu = np.sin((munu.mu - munu.node).to(u.radian).value)
    cosmu = np.cos((munu.mu - munu.node).to(u.radian).value)
    xx = cosmu * cosnu
    yy = sinmu * cosnu * cosi - sinnu * sini
    zz = sinmu * cosnu * sini + sinnu * cosi
    ra = Angle(np.arctan2(yy,xx),unit=u.radian) + munu.node
    dec = Angle(np.arcsin(zz),unit=u.radian)
    # if 'phi' in kwargs:
    #     phi = np.rad2deg(np.arctan2(cosmu * sini,
    #         (-sinmu * sinnu * sini + cosnu * cosi)*cosnu))
    #     return (ra,dec,phi)
    # else:
    #     return (ra,dec)
    return ICRS(ra=ra,dec=dec).transform_to(icrs_frame)
