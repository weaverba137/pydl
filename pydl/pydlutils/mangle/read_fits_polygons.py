#
#
#
def read_fits_polygons(filename,**kwargs):
    """Read a "polygon" format fits file.

    The main point of this is to extract the xcaps and cmcaps columns
    and replace them with caps.x and caps.cm, though this is really very silly
    & is only here so it's similar to IDL.
    """
    from pydlutils.goddard.fits import mrdfits
    poly,foo = mrdfits(filename,1,**kwargs)
    poly['CAPS'] = {'CM':list(), 'X':list()}
    for k in range(len(poly['NCAPS'])):
        poly['CAPS']['CM'].append(poly['CMCAPS'][k,0:poly['NCAPS'][k]])
        poly['CAPS']['X'].append(poly['XCAPS'][k].reshape((8,3))[0:poly['NCAPS'][k]])
    del poly['CMCAPS']
    del poly['XCAPS']
    return poly

