.. _pydl.photoop:

===================================
SDSS Imaging Data  (`pydl.photoop`)
===================================

Introduction
++++++++++++

The photoop_ package is used to process SDSS imaging data.  This package is used
to reduce_, resolve_ and flux-calibrate_ the SDSS raw data, resulting in both
flux-calibrated images and catalogs.

The primary emphasis of this implementation is:

1. Functions related to the SDSS photometric ``objID``, which is a unique
   integer used in SDSS databases, constructed from quantities that specify
   a particular astronomical object on a particular image.
2. Functions related to finding SDSS photometric data on disk.
3. The SDSS "window" function, which defines what parts of the sky are covered
   by SDSS images.  This can be used in conjunction with the
   :mod:`~pydl.pydlutils.mangle` module to find points and regions that have
   SDSS imaging.

The photoop package is itself divided into a number of subpackages.  Below we list
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
    Pretty much anything you could do with the photoop code you can do with the equivalent here.

=========== =============== ===================================================
Subpackage  Readiness Level Comments
=========== =============== ===================================================
apache      None
astrom      None
atlas       None
bluetip     None
compare     None
database    None
flats       None
hoggpipe    None
image       None
ircam       None
match       None
misc        None
pcalib      None
photoobj    Rudimentary
plan        None
plots       None
psf         None
ptcalib     None
resolve     None
sdss3_runqa None
sdssio      Rudimentary
window      Rudimentary
=========== =============== ===================================================

.. _photoop: https://svn.sdss.org/public/repo/sdss/photoop/trunk/
.. _reduce: https://www.sdss.org/dr14/imaging/pipeline/
.. _resolve: https://www.sdss.org/dr14/algorithms/resolve/
.. _flux-calibrate: https://www.sdss.org/dr14/algorithms/fluxcal/

API
+++

.. automodapi:: pydl.photoop

.. automodapi:: pydl.photoop.photoobj

.. automodapi:: pydl.photoop.sdssio

.. automodapi:: pydl.photoop.window
