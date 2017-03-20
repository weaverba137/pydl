.. _pydl.pydlutils:

=================================
SDSS Utilities (`pydl.pydlutils`)
=================================

Introduction
++++++++++++

This package provides functionality similar to idlutils_,
a general suite of tools heavily used by SDSS_.

.. _idlutils: http://www.sdss.org/dr12/software/idlutils/
.. _SDSS: http://www.sdss.org

idlutils_ is itself divided into a number of subpackages.  Below we list
the subpackages and the usability of the PyDL equivalent.
The readiness levels are defined as:

Not Applicable (NA)
    Another built-in or numpy/scipy/astropy package completely replaces this.  No point in implementing.
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
bspline     Good
cooling     Good
coord       Fair            Some functionality already provided by :mod:`astropy.coordinates`.
cosmography None
dimage      None
djsphot     None
dust        None            For use with the SFD galactic dust map.
first       None            For use with matching FIRST catalogs to SDSS data.
fits        NA              Use :mod:`astropy.io.fits`.
healpix     None
image       Rudimentary     Only a few functions implemented.
json        NA              Use :mod:`json` or other packages.
mangle      Fair            Some work still required on polygon area calculations.
math        Fair
mcmc        None            But there are plenty of good Python MCMC packages out there.
mglib       None
misc        Fair
mpeg        None
mpfit       None
physics     None
plot        None            Much functionality already exists in matplotlib.
psf         None
rgbcolor    Good
rosat       None            For use with matching ROSAT catalogs to SDSS data.
sdss        Good            Most important functionalities are bitmasks and reading sweep files.
slatec      None
spheregroup Good            Used for matching arbitrary RA, Dec coordinates to other arbitrary RA, Dec coordinates.
TeXtoIDL    NA              This package is for including TeX in IDL plots.  Since matplotlib understands TeX natively, this is not needed.
trace       Fair
ukidss      None            Used for matching UKIDSS catalogs to SDSS data.
wise        None            Used for matching WISE catalogs to SDSS data.
yanny       Good
=========== =============== ===================================================



API
+++

.. automodapi:: pydl.pydlutils

.. automodapi:: pydl.pydlutils.bspline

.. automodapi:: pydl.pydlutils.cooling

.. automodapi:: pydl.pydlutils.coord

.. automodapi:: pydl.pydlutils.image

.. automodapi:: pydl.pydlutils.mangle
    :skip: PydlutilsException, PydlutilsUserWarning

.. automodapi:: pydl.pydlutils.math

.. automodapi:: pydl.pydlutils.misc

.. automodapi:: pydl.pydlutils.rgbcolor

.. automodapi:: pydl.pydlutils.sdss

.. automodapi:: pydl.pydlutils.spheregroup
    :skip: PydlutilsException, PydlutilsUserWarning

.. automodapi:: pydl.pydlutils.trace
    :skip: PydlutilsException, flegendre

.. automodapi:: pydl.pydlutils.yanny
    :skip: OrderedDict, PydlutilsException, PydlutilsUserWarning, Table
