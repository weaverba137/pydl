# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the coord directory in idlutils.
"""
import numpy as np
from astropy import units as u
from astropy.coordinates import (Angle, BaseCoordinateFrame, FrameAttribute,
                                 SphericalRepresentation,
                                 QuantityFrameAttribute,
                                 RepresentationMapping,
                                 frame_transform_graph,
                                 FunctionTransform, ICRS)
from ..goddard.astro import get_juldate


class SDSSMuNu(BaseCoordinateFrame):
    """SDSS Great Circle Coordinates

    Attributes
    ----------
    stripe
        SDSS Stripe Number
    node
        Node of the great circle with respect to the celestial equator.
        In SDSS, this is almost always RA = 95.0 degrees.
    incl
        Inclination of the great circle with respect to the celestial
        equator.
    phi
        Counter-clockwise position angle w.r.t. north for an arc
        in the +nu direction.

    Parameters
    ----------
    mu : :class:`~astropy.coordinates.Angle`
        Angle corresponding to longitude measured along a stripe.
    nu : :class:`~astropy.coordinates.Angle`
        Angle corresponding to latitude measured perpendicular to a stripe.

    Notes
    -----
    http://www.sdss.org/dr12/algorithms/surveycoords/
    """
    default_representation = SphericalRepresentation
    frame_specific_representation_info = {
        'spherical': [
            RepresentationMapping(reprname='lon', framename='mu',
                                    defaultunit=u.deg),
            RepresentationMapping(reprname='lat', framename='nu',
                                    defaultunit=u.deg)
            ]
        }
    frame_specific_representation_info['unitspherical'] = (
            frame_specific_representation_info['spherical'])
    stripe = FrameAttribute(default=0)
    node = QuantityFrameAttribute(default=Angle(95.0, unit=u.deg), unit=u.deg)
    # phi = QuantityFrameAttribute(default=None, unit=u.deg)

    @property
    def incl(self):
        return Angle(stripe_to_incl(self.stripe), unit=u.deg)


def current_mjd():
    """Return the current MJD.
    """
    return get_juldate() - 2400000.5


@frame_transform_graph.transform(FunctionTransform, SDSSMuNu, ICRS)
def munu_to_radec(munu, icrs_frame):
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
    ra = Angle(np.arctan2(yy, xx), unit=u.radian) + munu.node
    dec = Angle(np.arcsin(zz), unit=u.radian)
    # if 'phi' in kwargs:
    #     phi = np.rad2deg(np.arctan2(cosmu * sini,
    #         (-sinmu * sinnu * sini + cosnu * cosi)*cosnu))
    #     return (ra, dec, phi)
    # else:
    #     return (ra, dec)
    return ICRS(ra=ra, dec=dec).transform_to(icrs_frame)


@frame_transform_graph.transform(FunctionTransform, ICRS, SDSSMuNu)
def radec_to_munu(icrs_frame, munu):
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
    mu = Angle(np.arctan2(y2, x2), unit=u.radian) + munu.node
    nu = Angle(np.arcsin(z2), unit=u.radian)
    # if 'phi' in kwargs:
    #     sinnu = np.sin(np.deg2rad(nu))
    #     cosnu = np.cos(np.deg2rad(nu))
    #     sinmu = np.sin(np.deg2rad(mu-node))
    #     cosmu = np.cos(np.deg2rad(mu-node))
    #     phi = np.rad2deg(np.arctan2(cosmu * sini,
    #         (-sinmu * sinnu * sini + cosnu * cosi)*cosnu))
    #     return (ra, dec, phi)
    # else:
    #     return (ra, dec)
    return SDSSMuNu(mu=mu, nu=nu, stripe=munu.stripe)


def stripe_to_eta(stripe):
    """Convert from SDSS great circle coordinates to equatorial coordinates.

    Parameters
    ----------
    stripe : :class:`int` or :class:`numpy.ndarray`
        SDSS Stripe number.

    Returns
    -------
    stripe_to_eta : :class:`float` or :class:`numpy.ndarray`
        The eta value in the SDSS (lambda,eta) coordinate system.
    """
    stripe_sep = 2.5
    eta = stripe * stripe_sep - 57.5
    if stripe > 46:
        eta -= 180.0
    return eta


def stripe_to_incl(stripe):
    """Convert from SDSS stripe number to an inclination relative to the
    equator.

    Parameters
    ----------
    stripe : :class:`int` or :class:`numpy.ndarray`
        SDSS Stripe number.

    Returns
    -------
    stripe_to_incl : :class:`float` or :class:`numpy.ndarray`
        Inclination of the stripe relative to the equator (Dec = 0).
    """
    dec_center = 32.5
    eta_center = stripe_to_eta(stripe)
    incl = eta_center + dec_center
    return incl
