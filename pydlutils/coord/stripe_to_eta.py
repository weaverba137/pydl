#
#
#
def stripe_to_eta(stripe):
    """Convert from SDSS great circle coordinates to equatorial coordinates.
    """
    stripe_sep = 2.5
    eta = stripe * stripe_sep - 57.5
    if stripe > 46:
        eta -= 180.0
    return eta

