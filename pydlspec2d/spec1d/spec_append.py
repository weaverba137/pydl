#
#
#
def spec_append(spec1, spec2, **kwargs):
    """Append the array spec2 to the array spec1 & return a new array.

    If the dimension of these arrays is the same, then append as [spec1,spec2].
    If not, increase the size of the smaller array & fill with zeros.
    If pixshift is set, the zeros will be prepended, otherwise they will be
    appended.
    """
    import numpy as np
    nrows1,npix1 = spec1.shape
    nrows2,npix2 = spec2.shape
    nrows = nrows1+nrows2
    nadd1 = 0
    nadd2 = 0
    if 'pixshift' in kwargs:
        if kwargs['pixshift'] < 0:
            nadd1 = -kwargs['pixshift']
        elif kwargs['pixshift'] > 0:
            nadd2 = kwargs['pixshift']
        else:
            pass
    maxpix = max(npix1+nadd1,npix2+nadd2)
    spec3 = np.zeros((nrows,maxpix),dtype=spec1.dtype)
    spec3[0:nrows1,nadd1:nadd1+npix1] = spec1
    spec3[nrows1:nrows,nadd2:nadd2+npix2] = spec2
    return spec3

