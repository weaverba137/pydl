#
#
#
def get_juldate():
    """Returns the current Julian date.

    Uses the MJD trick & adds the offset to get JD.
    """
    import time
    mjd = time.time()/86400.0 + 40587.0
    return mjd + 2400000.5

if __name__ == '__main__':
    print get_juldate()

