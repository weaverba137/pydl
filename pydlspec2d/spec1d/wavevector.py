#
#
#
def wavevector(minfullwave,maxfullwave,**kwargs):
    """Simply return an array of wavelengths.
    """
    import numpy as np
    if 'zeropoint' in kwargs:
        zeropoint = kwargs['zeropoint']
    else:
        zeropoint = 3.5
    if 'binsz' in kwargs:
        binsz = kwargs['binsz']
    else:
        binsz = 1.0e-4
    if 'wavemin' in kwargs:
        wavemin = kwargs['wavemin']
        spotmin = 0
        if 'wavemax' in kwargs:
            wavemax = kwargs['wavemax']
            spotmax = int((wavemax - wavemin)/binsz)
        else:
            spotmax = int((maxfullwave - wavemin)/binsz)
            wavemax = spotmax * binsz + wavemin
    else:
        spotmin = int((minfullwave - zeropoint)/binsz) + 1
        spotmax = int((maxfullwave - zeropoint)/binsz)
        wavemin = spotmin * binsz + zeropoint
        wavemax = spotmax * binsz + zeropoint
    nfinalpix = spotmax - spotmin + 1
    finalwave = np.arange(nfinalpix,dtype='d') *binsz + wavemin
    return finalwave

