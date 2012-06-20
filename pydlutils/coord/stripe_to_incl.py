#
#
#
def stripe_to_incl(stripe):
    """Convert from SDSS stripe number to an inclination relative to the equator.
    """
    from pydlutils.coord import stripe_to_eta
    dec_center = 32.5
    eta_center = stripe_to_eta(stripe)
    incl = eta_center + dec_center
    return incl

