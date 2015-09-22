# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def stripe_to_incl(stripe):
    """Convert from SDSS stripe number to an inclination relative to the equator.

    Parameters
    ----------
    stripe : :class:`int` or :class:`numpy.ndarray`
        SDSS Stripe number.

    Returns
    -------
    stripe_to_incl : :class:`float` or :class:`numpy.ndarray`
        Inclination of the stripe relative to the equator (Dec = 0).
    """
    from . import stripe_to_eta
    dec_center = 32.5
    eta_center = stripe_to_eta(stripe)
    incl = eta_center + dec_center
    return incl
