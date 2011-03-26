#
#
#
def hogg_mrdfits(filename,hdu=0,**kwargs):
    """Read a FITS file & return the data.

    Only bintables are supported right now.

    Note: pydlutils.goddard.fits.mrdfits is the 'real' port of
    hogg_mrdfits at least for now.  This function provides
    a placeholder for further development.
    """
    from pydlutils.goddard.fits import mrdfits
    return mrdfits(filename,hdu,**kwargs)
