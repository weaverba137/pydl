# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
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
