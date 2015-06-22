# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def airtovac(air):
    """Convert air wavelengths to wavelengths in vacuum.

    Parameters
    ----------
    air : array-like
        Values of wavelength in air in Angstroms.

    Returns
    -------
    airtovac : array-like
        Values of wavelength in vacuum in Angstroms.

    Notes
    -----
    * Formula from Ciddor 1996  Applied Optics , 35, 1566.
    * Values of wavelength below 2000 A are not converted.
    """
    from numpy import zeros
    try:
        vacuum = zeros(air.size,dtype=air.dtype) + air
        g = vacuum < 2000.0
    except AttributeError:
        # Most likely, vacuum is simply a float.
        vacuum = air
        g = None
        if air < 2000.0:
            return air
    for k in range(2):
        sigma2 = (1.0e4/vacuum)**2
        fact = 1.0 +  5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/( 57.362 - sigma2)
        vacuum = air * fact
    if g is not None:
        vacuum[g] = air[g]
    return vacuum
