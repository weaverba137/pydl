.. _pydl.photoop:

===================================
SDSS Imaging Data  (`pydl.photoop`)
===================================

Introduction
++++++++++++

The photoop_ package is used to process SDSS imaging data.  This package is used
to reduce_, resolve_ and flux-calibrate_ the SDSS raw data, resulting in both
flux-calibrated images and catalogs.

SDSS ceased taking imaging data in 2009, and there has only been one full
processing of the imaging data since then, although adjustments have been
made to the astrometry and flux calibration.

The primary emphasis of this implementation is:

1. Functions related to the SDSS photometric ``objID``, which is a unique
   integer used in SDSS databases, constructed from quantities that specify
   a particular astronomical object on a particular image.
2. Functions related to finding SDSS photometric data on disk.
3. The SDSS "window" function, which defines what parts of the sky are covered
   by SDSS images.  This can be used in conjunction with the
   :mod:`~pydl.pydlutils.mangle` module to find points and regions that have
   SDSS imaging.
4. In general, functions that work with existing imaging data, rather than
   functions to reduce the data.

The photoop package is itself divided into a number of subpackages.  Below we list
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
    Pretty much anything you could do with the photoop code you can do with the equivalent here.

=========== =============== ===================================================
Subpackage  Readiness Level Comments
=========== =============== ===================================================
apache      Obsolete        Processing of "Apache Wheel" images, used for calibration.
astrom      None            Astrometry for SDSS images.
atlas       None            Construction of "atlas" images, small cutouts of individual objects.
bluetip     None            Tools for "blue-tip" photometry and extinction estimation.
compare     Obsolete        Compare the same images in two different data reductions.
database    Obsolete        Experimental database loading code.
flats       None            Analysis of flat-field files.
hoggpipe    Obsolete        Another version of "Apache Wheel" processing code.
image       Rudimentary     Tools for creating "`corrected frame`_" images.  These are flux-calibrated and sky-subtracted images with physical flux units.
ircam       Obsolete        Tools for processing all-sky "cloud camera" images, used to establish photometricity.
match       None            Tools for matching SDSS spectra to corresponding photometric objects.
misc        None            Code with no obvious home in any other category.
pcalib      None            Tools related to "ubercalibration_".
photoobj    Rudimentary     Tools for creating calibrated catalogs from images.
plan        Obsolete        Tools for planning photometric reductions.
plots       None            Plots for high-level quality assurance.
psf         None            Analysis of point-spread functions.
ptcalib     Obsolete        Processing of "photometric telescope" data, an obsolete technique for flux-calibration.
resolve     None            Code for the resolve_ stage of image processing.
sdss3_runqa Obsolete        Quality assurance tests from the most recent photometric reduction.
sdssio      Rudimentary     Tools for reading and writing various data files produced by the photometric reductions.
window      Rudimentary     Tools for determining the sky coverage of the survey.
=========== =============== ===================================================

.. _photoop: https://svn.sdss.org/public/repo/sdss/photoop/trunk/
.. _reduce: https://www.sdss4.org/dr16/imaging/pipeline/
.. _resolve: https://www.sdss4.org/dr16/algorithms/resolve/
.. _flux-calibrate: https://www.sdss4.org/dr16/algorithms/fluxcal/
.. _`corrected frame`: https://data.sdss.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
.. _ubercalibration: https://www.sdss4.org/dr16/algorithms/fluxcal/

API
+++

.. automodapi:: pydl.photoop

.. automodapi:: pydl.photoop.image

.. automodapi:: pydl.photoop.photoobj

.. automodapi:: pydl.photoop.sdssio

.. automodapi:: pydl.photoop.window
    :skip: warn, sdss_name, sdss_calib, sdss_flagval, FITS_polygon
