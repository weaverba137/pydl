#
#
#
def is_cap_used(use_caps,i):
    """Returns whether a cap is used.
    """
    return (use_caps & 2L**i) > 0
