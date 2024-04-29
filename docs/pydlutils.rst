.. _pydl.pydlutils:

=================================
SDSS Utilities (`pydl.pydlutils`)
=================================

Introduction
++++++++++++

This package provides functionality similar to idlutils_,
a general suite of tools heavily used by SDSS_.

idlutils_ is itself divided into a number of subpackages.  Below we list
the subpackages and the usability of the PyDL equivalent.
The readiness levels are defined as:

Obsolete
    No point in implementing because the purpose of the code lapsed many years ago.
Not Applicable (NA)
    No point in implementing because another built-in or numpy/scipy/astropy package completely replaces this.
None
    Not (yet) implemented at all.
Rudimentary
    Only a few functions are implemented.
Fair
    Enough functions are implemented to be useful, but some are missing.
Good
    Pretty much anything you could do with the idlutils code you can do with the equivalent here.

=========== =============== ===================================================
Subpackage  Readiness Level Comments
=========== =============== ===================================================
2mass       None            For use with matching 2MASS catalogs to SDSS data.
astrom      None            For use with SDSS astrometric data structures.  Largely superseded by WCS.
bspline     Good            Fitting B-splines to data, especially for resampling.
cooling     Good            See :func:`pydl.pydlutils.cooling.read_ds_cooling`.
coord       Fair            Some functionality already provided by :mod:`astropy.coordinates`.
cosmography NA              Tools for computing lookback time, angular sizes at cosmological distances, etc. Use :mod:`astropy.cosmology`.
dimage      None            Interface to C code used for sky subtraction.
djsphot     None            A simple aperture photometry code.
dust        None            For use with the SFD galactic dust map.
first       None            For use with matching FIRST catalogs to SDSS data.
fits        NA              Use :mod:`astropy.io.fits`.
healpix     NA              Interact with HEALPix data.  Use healpy_.
image       Rudimentary     Image manipulation functions.
json        NA              Use :mod:`json` or other packages.
mangle      Fair            Some work still required on polygon area calculations.
math        Fair            Generic mathematical functions.  Many are implemented in numpy or scipy.
mcmc        None            But there are plenty of good Python MCMC packages out there.
mglib       Obsolete        An IDL object-oriented configuration file reader.
misc        Fair            General purpose utility functions.
mpeg        None            Wrapper for :command:`ppmtompeg`, makes movies from data.
mpfit       None            Appears to be an out-of-sync copy of the "markwardt" package in the `The IDL® Astronomy User's Libary`_.
physics     None            Implementation of physical formulas, *e.g.* free-free scattering.
plot        None            Much functionality already exists in matplotlib.
psf         Obsolete        Point-spread function fitting.
rgbcolor    Good            Some functionality is duplicated in :mod:`astropy.visualization`, especially :func:`~astropy.visualization.make_lupton_rgb`.
rosat       None            For use with matching ROSAT catalogs to SDSS data.
sdss        Good            Most important functionalities are bitmasks_ and reading `sweep files`_.
slatec      None            Fit B-splines using C code.
spheregroup Good            Used for matching arbitrary RA, Dec coordinates to other arbitrary RA, Dec coordinates.
TeXtoIDL    NA              This package is for including TeX in IDL plots.  Since matplotlib understands TeX natively, this is not needed.
trace       Fair            Used for fitting orthogonal functions to spectroscopic wavelength solutions.
ukidss      None            Used for matching UKIDSS catalogs to SDSS data.
wise        None            Used for matching WISE catalogs to SDSS data.
yanny       Good            Tools for manipulating `SDSS parameter files`_.
=========== =============== ===================================================

.. _idlutils: https://www.sdss4.org/dr16/software/idlutils/
.. _SDSS: https://www.sdss.org
.. _`The IDL® Astronomy User's Libary`: https://asd.gsfc.nasa.gov/archive/idlastro/
.. _healpy: https://healpy.readthedocs.io/en/latest/
.. _bitmasks: https://www.sdss4.org/dr16/algorithms/bitmasks/
.. _`sweep files`: https://data.sdss.org/datamodel/files/PHOTO_SWEEP/RERUN/calibObj.html
.. _`SDSS parameter files`: https://www.sdss4.org/dr16/software/par/

API
+++

.. automodapi:: pydl.pydlutils

.. automodapi:: pydl.pydlutils.bspline
    :skip: warn, PydlutilsUserWarning, djs_reject, fchebyshev, uniq, flegendre, cholesky_banded, LinAlgError, cho_solve_banded

.. automodapi:: pydl.pydlutils.cooling
    :skip: interp, get_pkg_data_contents

.. automodapi:: pydl.pydlutils.coord
    :skip: get_juldate

.. automodapi:: pydl.pydlutils.image

.. automodapi:: pydl.pydlutils.mangle
    :skip: PydlutilsException, PydlutilsUserWarning

.. automodapi:: pydl.pydlutils.math
    :skip: svd, djs_laxisnum, median

.. automodapi:: pydl.pydlutils.misc

.. automodapi:: pydl.pydlutils.rgbcolor
    :skip: warn

.. automodapi:: pydl.pydlutils.sdss
    :skip: download_file, spherematch, uniq

.. automodapi:: pydl.pydlutils.spheregroup
    :skip: warn, PydlutilsException, PydlutilsUserWarning, gcirc

.. automodapi:: pydl.pydlutils.trace
    :skip: chebyt, FITS_rec, PydlutilsException, djs_reject, djs_laxisgen, flegendre

.. automodapi:: pydl.pydlutils.yanny
    :skip: OrderedDict, PydlutilsException, PydlutilsUserWarning, Table
