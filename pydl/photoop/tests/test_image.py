# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test the functions in pydl.photoop.image.
"""
import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from ..image import sdss_psf_recon


def test_sdss_psf_recon():
    psField = get_pkg_data_filename('t/psField-003366-3-0110.fit')
    with fits.open(psField) as hdulist:
        psf = sdss_psf_recon(hdulist[3].data, 600, 500)
    assert np.allclose(psf[23:26, 25],
                       np.array([0.02220068,  0.0738798 ,  0.11940149]))


def test_sdss_psf_norm():
    psField = get_pkg_data_filename('t/psField-003366-3-0110.fit')
    with fits.open(psField) as hdulist:
        psf = sdss_psf_recon(hdulist[3].data, 600, 500, normalize=100.0)
    assert np.allclose(psf.sum(), 100.0)
    assert np.allclose(psf[23:26, 25],
                       np.array([2.194198,   7.301887,  11.801009]))


def test_sdss_psf_recon_trim():
    psField = get_pkg_data_filename('t/psField-003366-3-0110.fit')
    with fits.open(psField) as hdulist:
        psf = sdss_psf_recon(hdulist[3].data, 600, 500, trimdim=(25, 25))
    # print(psf)
    assert np.allclose(psf[10:13, 12],
                       np.array([0.02220068,  0.0738798 ,  0.11940149]))
