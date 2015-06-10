# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
from astropy import units as u
from astropy.coordinates import (BaseCoordinateFrame, FrameAttribute,
    SphericalRepresentation, TimeFrameAttribute, QuantityFrameAttribute,
    RepresentationMapping)
#
class SDSSMuNuFrame(BaseCoordinateFrame):
    default_representation = SphericalRepresentation
    frame_specific_representation_info = {
        'spherical': [
            RepresentationMapping(reprname='lon', framename='mu', defaultunit=u.deg),
            RepresentationMapping(reprname='lat', framename='nu', defaultunit=u.deg),
            RepresentationMapping(reprname='distance', framename='DIST', defaultunit=None)
            ],
        'unitspherical': [
            RepresentationMapping(reprname='lon', framename='mu', defaultunit=u.deg),
            RepresentationMapping(reprname='lat', framename='nu', defaultunit=u.deg)
            ],
        'cartesian': [
            RepresentationMapping(reprname='x', framename='X'),
            RepresentationMapping(reprname='y', framename='Y'),
            RepresentationMapping(reprname='z', framename='Z')
            ]
        }
    stripe = FrameAttribute(default=None)
    node = QuantityFrameAttribute(default=None,unit=u.deg)
    incl = QuantityFrameAttribute(default=None,unit=u.deg)
    phi = QuantityFrameAttribute(default=None,unit=u.deg)
    equinox = TimeFrameAttribute(default='J2000')
