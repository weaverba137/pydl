# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the image directory in photoop.
"""
import numpy as np


def sdss_psf_recon(pstruct, row, col, counts=None, trim=False):
    """Reconstruct the PSF at position (`row`,`col`) in an SDSS field from the
    KL-decompositions.

    These can be read from the `psField files`_.

    .. _`psField files`: https://data.sdss.org/datamodel/files/PHOTO_REDUX/RERUN/RUN/objcs/CAMCOL/psField.html

    Parameters
    ----------
    pstruct
        A PSF KL-decomposition structure read from a psField file.
    row : :class:`int`
        Row in the field.
    col : :class:`int`
        Column in the field.
    counts : :class:`int` or :class:`float`, optional
        The total counts in the image. Default is whatever comes from the reconstruction, which is usually close to unity.
    trim : :class:`bool`, optional
        If ``True``, remove regions that contain only zeroes.

    Returns
    -------
    :class:`numpy.ndarray`
        The 2D reconstructed PSF.

    Example
    -------
    >>> import numpy as np
    >>> from astropy.io import fits
    >>> from astropy.utils.data import get_readable_fileobj
    >>> from pydl.photoop.image import sdss_psf_recon
    >>> psfile = ('https://data.sdss.org/sas/dr14/eboss/photo/redux/301/' +
    ...           '3366/objcs/3/psField-003366-3-0110.fit')
    >>> with get_readable_fileobj(psfile, encoding='binary') as psField:  # doctest: +REMOTE_DATA
    ...     with fits.open(psField) as hdulist:
    ...         psf = sdss_psf_recon(hdulist[3].data, 500, 500)
    """
    rcs = 0.001
    nrow_b = pstruct['nrow_b'][0]
    ncol_b = pstruct['ncol_b'][0]
    rnrow = pstruct['RNROW'][0]
    rncol = pstruct['RNCOL'][0]
    nb = nrow_b * ncol_b
    coeffs = np.zeros((nb,), dtype=np.float32)
    ecoeff = np.zeros((3,), dtype=np.float32)
    for i in range(nb):
        coeffs[i] = (row*rcs)**(i % nrow_b) * (col*rcs)**(i//nrow_b)
    for j in range(3):
        for i in range(nb):
            ecoeff[j] += coeffs[i] * pstruct['c'][j, i % nrow_b, i//nrow_b]
    p = pstruct['RROWS'][0]*ecoeff[0] + pstruct['RROWS'][1]*ecoeff[1] + pstruct['RROWS'][2]*ecoeff[2]
    if counts is not None:
        p /= p.sum()
        p *= counts
    p = p.reshape(rnrow, rncol)
    if trim:
        p = p[10:41, 10:41]
    return p
