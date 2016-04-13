# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the rgbcolor directory in idlutils.

This code reproduces the algorithms of Lupton *et al.* (2004) [2]_.
The purpose is to produce nice color (RGB, JPEG, PNG, etc.) from FITS
images in astronomical filter bands.

References
----------

..  [2] `Lupton, Robert, et al., 2004 PASP 116, 113
    <http://adsabs.harvard.edu/abs/2004PASP..116..133L>`_.
"""
import numpy as np


def nw_arcsinh(colors, nonlinearity=3.0):
    """Scales `colors` by a degree of nonlinearity specified by user.

    The input image must have zero background (*i.e.*, it must already
    be background-subtracted).

    Parameters
    ----------
    colors : :class:`~numpy.ndarray`
        3D Array containing RGB image.  The dimensions should be (X, Y, 3).
    nonlinearity : :class:`float`
        Amount of nonlinear scaling.  If set to zero, no scaling will be
        performed (this is equivalent to linear scaling).

    Returns
    -------
    :class:`~numpy.ndarray`
        The scaled image.

    Raises
    ------
    ValueError
        If `colors` has the wrong shape.
    """
    if nonlinearity == 0:
        return colors
    if colors.ndim != 3:
        raise ValueError("A 3D image is required!")
    nx, ny, nc = colors.shape
    if nc != 3:
        raise ValueError("The 3D image must have 3 image planes.")
    radius = colors.sum(2)
    w = radius == 0
    radius[w] = radius[w] + 1.0
    fac = np.arcsinh(radius*nonlinearity)/nonlinearity/radius
    fitted_colors = np.zeros(colors.shape, dtype=colors.dtype)
    for k in range(3):
        fitted_colors[:, :, k] = colors[:, :, k]*fac
    return fitted_colors


def nw_cut_to_box(colors, origin=(0.0, 0.0, 0.0)):
    """Limits the pixel values of the image to a 'box', so that the colors
    do not saturate to white but to a specific color.

    Parameters
    ----------
    colors : :class:`~numpy.ndarray`
        3D Array containing RGB image.  The dimensions should be (X, Y, 3).
    origin : :func:`tuple` or :class:`~numpy.ndarray`
        An array with 3 elements.  The "distance" from this origin is
        considered saturated.

    Returns
    -------
    :class:`~numpy.ndarray`
        The "boxed" image.

    Raises
    ------
    ValueError
        If `colors` or `origin` has the wrong shape.
    """
    if len(origin) != 3:
        raise ValueError("The origin array must contain 3 elements!")
    if colors.ndim != 3:
        raise ValueError("A 3D image is required!")
    nx, ny, nc = colors.shape
    if nc != 3:
        raise ValueError("The 3D image must have 3 image planes.")
    pos_dist = 1.0 - np.array(origin, dtype=colors.dtype)
    factors = np.zeros(colors.shape, dtype=colors.dtype)
    for k in range(3):
        factors[:, :, k] = colors[:, :, k]/pos_dist[k]
    factor = factors.max(2)
    factor = np.where(factor > 1.0, factor, 1.0)
    boxed_colors = np.zeros(colors.shape, dtype=colors.dtype)
    for k in range(3):
        boxed_colors[:, :, k] = origin[k] + colors[:, :, k]/factor
    return boxed_colors


def nw_float_to_byte(image, bits=8):
    """Converts an array of floats in [0.0, 1.0] to bytes in [0, 255].

    Parameters
    ----------
    image : :class:`~numpy.ndarray`
        Image to convert.
    bits : :class:`int`, optional
        Number of bits in final image.

    Returns
    -------
    :class:`~numpy.ndarray`
        Converted image.
    """
    from warnings import warn
    from . import PydlutilsUserWarning
    if bits > 8:
        warn("bits > 8 not supported; setting bits = 8.", PydlutilsUserWarning)
        bits = 8
    fmax = 1 << bits
    bmax = fmax - 1
    f1 = np.floor(image * fmax)
    f2 = np.where(f1 > 0, f1, 0)
    f3 = np.where(f2 < bmax, f2, bmax)
    return f3.astype(np.uint8)


def nw_scale_rgb(colors, scales=(1.0, 1.0, 1.0)):
    """Multiply RGB image by color-dependent scale factor.

    Parameters
    ----------
    colors : :class:`~numpy.ndarray`
        3D Array containing RGB image.  The dimensions should be (X, Y, 3).
    scales : :func:`tuple` or :class:`~numpy.ndarray`, optional
        An array with 3 elements.

    Returns
    -------
    :class:`~numpy.ndarray`
        The scaled image.

    Raises
    ------
    ValueError
        If `colors` or `scales` has the wrong shape.
    """
    if len(scales) != 3:
        raise ValueError("The scale factor must contain 3 elements!")
    if colors.ndim != 3:
        raise ValueError("A 3D image is required!")
    nx, ny, nc = colors.shape
    if nc != 3:
        raise ValueError("The 3D image must have 3 image planes.")
    scaled_colors = np.zeros(colors.shape, dtype=colors.dtype)
    for k in range(3):
        scaled_colors[:, :, k] = colors[:, :, k]*scales[k]
    return scaled_colors
