# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
def sdssflux2ab(flux,magnitude=False,ivar=False):
    """Convert the SDSS calibrated fluxes (magnitudes) into AB fluxes (magnitudes)

    Parameters
    ----------
    flux : ndarray
        Array of calibrated fluxes or SDSS magnitudes with 5 columns,
        corresponding to the 5 filters u,g,r,i,z
    magnitude : bool, optional
        If set to ``True``, then assume `flux` are SDSS magnitudes instead of linear
        flux units
    ivar : ndarray, optional
        If set, the input fluxes are actually inverse variances.

    Returns
    -------
    sdssflux2ab : ndarray
        Array of fluxes or magnitudes on the AB system.

    Notes
    -----
    Uses the conversions posted by D.Hogg (sdss-calib/845)::

     u(AB,2.5m) = u(2.5m) - 0.042
     g(AB,2.5m) = g(2.5m) + 0.036
     r(AB,2.5m) = r(2.5m) + 0.015
     i(AB,2.5m) = i(2.5m) + 0.013
     z(AB,2.5m) = z(2.5m) - 0.002
    """
    import numpy as np
    #
    # Correction vector, adjust this as necessary
    #
    correction = np.array([-0.042, 0.036, 0.015, 0.013, -0.002])
    rows, cols = flux.shape
    abflux = flux.copy()
    if magnitude:
        for i in range(rows):
            abflux[i,:] += correction
    else:
        factor = 10.0**(-correction/2.5)
        if ivar:
            factor = 1.0/factor**2
        for i in range(rows):
            abflux[i,:] *= factor
    return abflux
