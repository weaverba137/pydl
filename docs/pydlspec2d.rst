.. _pydl.pydlspec2d:

============================================
SDSS Spectroscopic Data  (`pydl.pydlspec2d`)
============================================

Introduction
++++++++++++

This package provides functionality in the SDSS idlspec2d_ package.  This
package is used for processing and analyzing data from the SDSS optical
spectrographs.  The code is thus relevant to the `SDSS Legacy`_,
BOSS_ and eBOSS_ surveys.  This package does *not* work with any infrared
spectrograph data associated with the APOGEE-2_ survey.

The primary *technical* focus of this particular implementation is the function
:func:`~pydl.pydlspec2d.spec2d.combine1fiber`.  This function is responsible
for resampling 1D spectra onto a new wavelength solution.  This allows for:

1. Shifting a spectrum from observed redshift to rest frame.
2. Coaddition of spectra of the same object, after resampling all spectra onto
   the same wavelength solution.

The primary *scientific* motivation of implementing
:func:`~pydl.pydlspec2d.spec2d.combine1fiber` is to create template spectra
based on curated spectra of, *e.g.*, luminous red galaxies (LRGs).
Principal Component Analysis (PCA) or other techniques may be used to
construct template spectra, but putting all spectra on the same
rest-frame wavelength solution is the first step.

The idlspec2d package is itself divided into a number of subpackages.  Below we list
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
    Pretty much anything you could do with the idlspec2d code you can do with the equivalent here.

========== =============== ===================================================
Subpackage Readiness Level Comments
========== =============== ===================================================
apo2d      None            Quick extraction code for quality assurance at observation time.
config     None            Extraction pipeline configuration parameters in object-oriented IDL.
fluxfix    None            Flux calibration
guider     None            Interface to guider camera.
inspect    None            Tools for manual inspection of spectra.
photoz     Obsolete        Photometric redshifts for SDSS objects using spectral templates.
plan       None            Tools for planning exposures and recordings summaries of exposures.
plate      Obsolete        Tools for designing SDSS spectroscopic plates, especially for star clusters.
science    None            Code for science analysis of sets of 1D spectra.
spec1d     Fair            Tools for processing 1D spectra, including redshift fitting.
spec2d     Fair            Tools for extracting spectra from 2D images.
specdb     Obsolete        Tests on storing spectroscopic results in SQL databases.
specflat   None            Flat-fielding of spectroscopic 2D images.
templates  None            Tools for constructing spectroscopic templates.
testsuite  None            Tools for high-level quality assurance, *e.g.* comparing two reductions of the same data.
========== =============== ===================================================

.. _idlspec2d: https://svn.sdss.org/public/repo/eboss/idlspec2d/trunk/
.. _`SDSS Legacy`: https://classic.sdss.org/legacy/index.html
.. _BOSS: https://www.sdss.org/surveys/boss/
.. _eBOSS: https://www.sdss.org/surveys/eboss/
.. _APOGEE-2: https://www.sdss.org/surveys/apogee-2/

API
+++

.. automodapi:: pydl.pydlspec2d

.. automodapi:: pydl.pydlspec2d.spec1d

.. automodapi:: pydl.pydlspec2d.spec2d
    :skip: warn, erf, get_pkg_data_filename, smooth, iterfit, djs_maskinterp, djs_median, sdss_flagval, traceset2xy, xy2traceset, vactoair
