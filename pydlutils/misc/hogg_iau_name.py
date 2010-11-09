#
#
#
def hogg_iau_name(ra,dec,prefix='SDSS',precision=1):
    """Properly format astronomical source names to the IAU convention.

    Arguments:
    ra -- Right ascencion in decimal degrees
    dec -- Declination in decimal degrees.
    """
    import numpy as np
    #
    # Promote scalar values to arrays.
    #
    if isinstance(ra,float):
        ra = np.array([ra])
    if isinstance(dec,float):
        dec = np.array([dec])
    h = ra/15.0
    rah = np.floor(h)
    ram = np.floor(60.0*(h-rah))
    ras = 60.0*(60.0*(h-rah) - ram)
    ras = np.floor(ras*10.0**(precision+1))/10.0**(precision+1)
    rasformat = "%%0%d.%df" % (precision+4, precision+1)
    desgn = np.array(list('+'*len(dec)))
    desgn[pylab.find(dec < 0)] = '-'
    adec = np.absolute(dec)
    ded = np.floor(adec)
    dem = np.floor(60.0*(adec-ded))
    des = 60.0*(60.0*(adec-ded) - dem)
    des = np.floor(des*10.0**precision)/10.0**precision
    desformat = "%%0%d.%df" % (precision+3, precision)
    if precision == 0:
        desformat = "%02d"
    adformat = "%%02d%%02d%s%%s%%02d%%02d%s" % (rasformat,desformat)
    adstr = map(lambda a,b,c,d,e,f,g: adformat % (a,b,c,d,e,f,g),
        rah, ram, ras, desgn, ded, dem, des)
    if prefix == '':
        jstr = 'J'
    else:
        jstr = ' J'
    name = map(lambda x: "%s%s%s" % (prefix, jstr, x), adstr)
    if len(ra) == 1:
        return name[0]
    else:
        return name

