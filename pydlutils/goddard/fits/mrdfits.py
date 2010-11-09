#
#
#
def mrdfits(filename,hdu=0,**kwargs):
    """Read a FITS file & return the data.

    Only bintables are supported right now.
    """
    import pyfits
    if 'range' in kwargs and 'rows' in kwargs:
        print "Keywords 'range' and 'rows' are mutually exclusive."
        return None
    try:
        fits = pyfits.open(filename)
    except IOError:
        try:
            fits = pyfits.open(filename+'.gz')
        except IOError:
            return None
    fits_cols = fits[hdu].columns.names
    fits_hdr = fits[hdu].header.ascardlist()
    data = dict()
    if 'columns' in kwargs:
        cols = kwargs['columns']
    else:
        cols = fits_cols
    if 'rows' in kwargs:
        rowdata = fits[hdu].data[kwargs['rows']]
        for col in cols:
            if col not in fits_cols:
                print "Invalid column: %s!" % col
            else:
                data[col] = rowdata[col]
    else:
        if 'range' in kwargs:
            if kwargs['range'] == 'all':
                r = (0,0)
            else:
                if isinstance(kwargs['range'],int):
                    r = (0,kwargs['range'])
                else:
                    r = kwargs['range']
        else:
            r = (0,0)
        for col in cols:
            if col not in fits_cols:
                print "Invalid column: %s!" % col
            else:
                if r[1] > 0:
                    data[col] = fits[hdu].data.field(col)[r[0]:r[1]]
                else:
                    data[col] = fits[hdu].data.field(col)
    fits.close()
    return (data,fits_hdr)
