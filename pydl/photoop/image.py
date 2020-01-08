# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the image directory in photoop.
"""
import numpy as np


def sdss_psf_recon(psfield, xpos, ypos, normalize=None, trimdim=None):
    """Reconstruct the PSF at position (`xpos`, `ypos`) in an SDSS field.

    These can be read from the `psField files`_.

    .. _`psField files`: https://data.sdss.org/datamodel/files/PHOTO_REDUX/RERUN/RUN/objcs/CAMCOL/psField.html

    Parameters
    ----------
    psfield : :class:`~astropy.io.fits.fitsrec.FITS_rec`
        A PSF KL-decomposition structure read from a psField file.  This is
        from HDU's 1 through 5 in a psField file, corresponding
        to the five filters *u*, *g*, *r*, *i*, *z*.
    xpos : :class:`int`
        Column position (0-indexed, not 0.5-indexed as PHOTO outputs).
    ypos : :class:`int`
        Row position (0-indexed, not 0.5-indexed as PHOTO outputs).
    normalize : :class:`int` or :class:`float`, optional
        If set, normalize the integral of the image to this value.
    trimdim : :class:`tuple`, optional
        Trimmed dimensions; for example, set to (25, 25) to trim
        the output PSF image to those dimensions.  These dimensions
        must be odd-valued.

    Returns
    -------
    :class:`numpy.ndarray`
        The 2D reconstructed PSF image, typically dimensioned (51, 51).
        The center of the PSF is always the central pixel; this function will
        not apply any sub-pixel shifts.

    Notes
    -----
    The SDSS photo PSF is described as a set of eigen-templates, where the
    mix of these eigen-templates is a simple polynomial function with (x,y)
    position on the CCD.  Typically, there are 4 such 51x51 pixel templates.
    The polynomial functions are typically quadratic in each dimension,
    with no cross-terms.

    The formula is the following, where :math:`i` is the index of row
    polynomial order, :math:`j` is the index of column polynomial order,
    and :math:`k` is the template index:

    .. math::

        a_k &= \\sum_i \\sum_j (0.001 \\times \\mathrm{ROWC})^i \\times (0.001 \\times \\mathrm{COLC})^j \\times [C_k]_{ij} \\\\
        \\mathrm{psfimage} &= \\sum_k a_k \\mathrm{RROWS}_k

    The polynomial terms need not be of the same order for each template.

    Examples
    --------
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
    #
    # Hard-wired scale factor for ypos/xpos coefficients.
    #
    rc_scale = 0.001
    #
    # Assume that the dimensions of each eigen-template are the same.
    #
    rnrow = psfield['RNROW'][0]
    rncol = psfield['RNCOL'][0]
    npix = rnrow * rncol
    #
    # These are the polynomial coefficients as a function of x,y.
    # Only compute these coefficients for the maximum polynomial order in use.
    # In general, this order can be different for each eigen-template.
    #
    nr_max = psfield['nrow_b'].max()
    nc_max = psfield['ncol_b'].max()
    # nb = nrow_b * ncol_b
    coeffs = np.outer(((ypos+0.5) * rc_scale)**np.arange(nr_max), ((xpos+0.5) * rc_scale)**np.arange(nc_max))
    #
    # Reconstruct the image by summing each eigen-template.
    #
    neigen = psfield.shape[0]
    psfimage = np.zeros(psfield['RROWS'][0].shape,
                        dtype=psfield['RROWS'][0].dtype)
    for i in range(neigen):
        psfimage += (psfield['RROWS'][i] *
                     (psfield[i]['c'][0:psfield[0]['nrow_b'], 0:psfield[0]['ncol_b']] *
                      coeffs[0:psfield[0]['nrow_b'], 0:psfield[0]['ncol_b']].T).sum())
    #
    # We have reconstructed the PSF as a vector using all the
    # pixels in psfield.RROWS. So, at the end we trim this vector to only
    # those pixels used, and reform it into a 2-dimensional image.
    #
    if normalize is not None:
        psfimage /= psfimage.sum()
        psfimage *= normalize
    psfimage = psfimage[0:npix].reshape(rncol, rnrow)
    if trimdim is not None:
        psfimage = psfimage[(rncol - trimdim[0])//2:(rncol + trimdim[0])//2,
                            (rnrow - trimdim[1])//2:(rnrow + trimdim[1])//2]
    return psfimage
