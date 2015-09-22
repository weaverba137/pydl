# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from astropy import units as u
from astropy.coordinates import (Angle, BaseCoordinateFrame, FrameAttribute,
    SphericalRepresentation, QuantityFrameAttribute,
    RepresentationMapping)
from .stripe_to_incl import stripe_to_incl
#
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
            RepresentationMapping(reprname='lon', framename='mu', defaultunit=u.deg),
            RepresentationMapping(reprname='lat', framename='nu', defaultunit=u.deg)
            ]
        }
    frame_specific_representation_info['unitspherical'] = frame_specific_representation_info['spherical']
    stripe = FrameAttribute(default=0)
    node = QuantityFrameAttribute(default=Angle(95.0,unit=u.deg),unit=u.deg)
    # phi = QuantityFrameAttribute(default=None,unit=u.deg)
    #
    #
    #
    @property
    def incl(self):
        return Angle(stripe_to_incl(self.stripe),unit=u.deg)
