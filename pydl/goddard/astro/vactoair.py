# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def vactoair(vacuum):
    """Convert vacuum wavelengths to wavelengths in air.

    Parameters
    ----------
    vacuum : array-like
        Values of wavelength in vacuum in Angstroms.

    Returns
    -------
    vactoair : array-like
        Values of wavelength in air in Angstroms.

    Notes
    -----
    * Formula from `P. E. Ciddor, Applied Optics, 35, 1566 (1996) <http://adsabs.harvard.edu/abs/1996ApOpt..35.1566C>`_.
    * Values of wavelength below 2000 Å are not converted.
    """
    from numpy import zeros
    try:
        air = zeros(vacuum.size,dtype=vacuum.dtype)
        g = vacuum < 2000.0
    except AttributeError:
        # Most likely, vacuum is simply a float.
        air = vacuum
        g = None
        if vacuum < 2000.0:
            return air
    sigma2 = (1.0e4/vacuum)**2
    fact = 1.0 +  5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/( 57.362 - sigma2)
    air = vacuum/fact
    if g is not None:
        air[g] = vacuum[g]
    return air
