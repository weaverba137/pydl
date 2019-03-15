==============
PyDL Changelog
==============

1.0.0 (unreleased)
------------------

*This version will only support Python 3 and Astropy 3.*  This release is
planned for early 2019.

* Remove all Python 2 constructs; add tests in :mod:`~pydl.photoop.window` (PR `#51`_).

.. _`#51`: https://github.com/weaverba137/pydl/pull/51

0.7.0 (2019-02-22)
------------------

*This version is planned to be the last version with Python 2 support.*

* Support the ``firstField`` bit in ObjIDs from DR7 and earlier (Issue `#37`_).
* Change tests of Astropy development version from Python 2 to Python 3.
* Update to `astropy_helpers`_/v2.0.6 (PR `#40`_).
* Add :mod:`astropy.units` support to :func:`~pydl.goddard.astro.airtovac`
  and :func:`~pydl.goddard.astro.vactoair` (PR `#41`_).
* Change Exelis to Harris Geospatial (PR `#42`_).
* Fix ``FutureWarning`` in ``re`` in Python 3.7 due to nested sets (PR `#44`_).
* Use ``six`` instead of ``astropy.extern.six`` (PR `#48`_).
* Update ``_astropy_init.py`` (PR `#47`_ via PR `#48`_).
* Update `astropy_helpers`_/v2.0.8 (PR `#45`_ via PR `#48`_).

.. _`#37`: https://github.com/weaverba137/pydl/issues/37.
.. _`#40`: https://github.com/weaverba137/pydl/pull/40
.. _`#41`: https://github.com/weaverba137/pydl/pull/41
.. _`#42`: https://github.com/weaverba137/pydl/pull/42
.. _`#44`: https://github.com/weaverba137/pydl/pull/44
.. _`#45`: https://github.com/weaverba137/pydl/pull/45
.. _`#47`: https://github.com/weaverba137/pydl/pull/47
.. _`#48`: https://github.com/weaverba137/pydl/pull/48

0.6.0 (2017-09-19)
------------------

* This release is compatible with Astropy 2.0, and may be backwards
  incompatible with astropy v1.x.
* Update to `astropy_helpers`_/v2.0.1.
* Use standard library :mod:`argparse` (Issue `#31`_).
* Use the new :class:`astropy.coordinates.Attribute` class.
* Fix typo (PR `#26`_).

.. _`#31`: https://github.com/weaverba137/pydl/issues/31.
.. _`#26`: https://github.com/weaverba137/pydl/pull/26

0.5.4 (2017-05-04)
------------------

* Added :func:`~pydl.pydlutils.sdss.sdss_specobjid` to compute SDSS
  specObjIDs, and its inverse function
  :func:`~pydl.pydlutils.sdss.unwrap_specobjid`.
* Update to `astropy_helpers`_/v1.3.1.
* Refactor HMF code into an object to contain the data and methods.
* Use functions from :mod:`astropy.utils.data` where possible.
* Fix an integer division error encountered when using Numpy 1.12
  (Issue `#19`_).
* Fixed tests that were failing on 32-bit platforms *and* Python 3.5
  (Issue `#20`_).

.. _`#19`: https://github.com/weaverba137/pydl/issues/19
.. _`#20`: https://github.com/weaverba137/pydl/issues/20

0.5.3 (2016-12-03)
------------------

* Fixed formatting of TODO document.
* Fixed tests that were failing on 32-bit platforms (Issue `#14`_).
* Use temporary files so that tests can run when astropy is installed
  read-only (*e.g.*, with :command:`pip`; Issue `#16`_)

.. _`#14`: https://github.com/weaverba137/pydl/issues/14
.. _`#16`: https://github.com/weaverba137/pydl/issues/16

0.5.2 (2016-08-04)
------------------

* Changes in how Mangle-polygon containing FITS files are handled, related to
  Issue `#11`_.
* Fixed memory leak in :func:`~pydl.pydlspec2d.spec2d.combine1fiber`,
  see Issue `#12`_.
* Added :func:`~pydl.pydlutils.mangle.is_in_window`.
* Allow polygon area functions to deal with negative caps and ``use_caps``.
* Update ``docs/conf.py`` for Python 3.5 compatibility (PR `#13`_).

.. _`#13`: https://github.com/weaverba137/pydl/pull/13
.. _`#11`: https://github.com/weaverba137/pydl/issues/11
.. _`#12`: https://github.com/weaverba137/pydl/issues/12


0.5.1 (2016-06-22)
------------------

* Removed unnecessary ``from __future__`` import in
  :mod:`pydl.pydlspec2d.spec1d`.
* Ongoing documentation upgrades.
* Update some links that needed to be transitioned from SDSS-III to SDSS-IV.
* Upgrade to `astropy_helpers`_/v1.2.
* Update to latest version of package-template_.
* Disabled tests on Python 3.3; enabled tests on Python 3.5
* Fix Issue `#8`_; Issue `#9`_.
* Add warnings about incomplete Mangle functions.

.. _`#8`: https://github.com/weaverba137/pydl/issues/8
.. _`#9`: https://github.com/weaverba137/pydl/issues/9

0.5.0 (2016-05-01)
------------------

* Dropped support for Python 2.6.  Python 2.6 does not contain
  :class:`collections.OrderedDict`, which is needed to support
  :class:`~pydl.pydlutils.yanny.yanny` objects, and at this point it is not
  worth going to the trouble to support this with an external package.
* Ongoing review and upgrade of docstrings.
* Yanny files can now be converted into *genuine* NumPy
  :class:`record arrays <numpy.recarray>`; previously, the conversion was only
  to :class:`numpy.ndarray` with named columns, which is a slightly different
  thing.
* Added additional tests on :class:`~pydl.pydlutils.yanny.yanny` objects.
* Experimental support for interconversion of
  :class:`~pydl.pydlutils.yanny.yanny` objects and
  :class:`~astropy.table.Table` objects.
* Improving `PEP 8`_ compliance
* Restructuing sub-packages to reduce the number of files.
* Improvements to spectral template processing code, deduplicated some code.
* Support platform-independent home directory (PR `#7`_).
* Uppercase the package name (in documentation only).
* Upgrade to `astropy_helpers`_/v1.1.1.
* Add functions from the idlutils rgbcolor directory.
* :func:`~pydl.pydlspec2d.spec1d.spec_path` can now find SDSS spectra, not just
  BOSS.

.. _`PEP 8`: https://www.python.org/dev/peps/pep-0008/
.. _`#7`: https://github.com/weaverba137/pydl/pull/7

0.4.1 (2015-09-22)
------------------

* No changes at all from 0.4.0.  This tag only exists because of a botched
  PyPI upload.

0.4.0 (2015-09-22)
------------------

* Use `astropy_helpers`_/v1.0.3, package-template_/v1.0.
* Remove some old FITS code that :mod:`astropy.io.fits` makes moot.
* Remove code for command-line scripts.  These are now auto-generated by the
  "entry_point" method.
* Remove Python/3.2 tests.
* Improved test coverage.
* Fixed problem with the :mod:`~pydl.pydlutils.spheregroup` code.
* Removed some code that is 100% redundant with astropy (*e.g.* ``angles_to_xyz()``).
* Fixed bug in :func:`~pydl.pydlutils.mangle.set_use_caps` that was discovered on the IDL side.
* Updated documentation of :func:`~pydl.pydlutils.mangle.read_fits_polygons`.
* Added cross-references to classes, functions, etc.

0.3.0 (2015-02-20)
------------------

* Use `astropy_helpers`_/v0.4.3, package-template_/v0.4.1.
* Avoided (but did not fix) a bug in :class:`~pydl.pydlutils.spheregroup.chunks` that occurs when operating on
  a list of coordinates of length 1.
* Fixed a typo in :class:`~pydl.pydlutils.bspline.bspline`, added documentation.
* Simplify documentation files.
* :func:`~pydl.pydlutils.sdss.sdss_flagname` now accepts more types of numeric input.
* Added :doc:`credits` file.

0.2.3 (2014-07-22)
------------------

* Added :mod:`pydl.photoop.window`.
* Added stub :func:`~pydl.photoop.sdssio.sdss_calib`, updated :func:`~pydl.photoop.window.sdss_score`.
* Added :func:`~pydl.photoop.photoobj.unwrap_objid`.
* Merged pull request #4, fixing some Python3 issues.

0.2.2 (2014-05-07)
------------------

* Updated to latest package-template_ version.
* Added ability to `write multiple ndarray to yanny files`_.
* Fixed :func:`~pydl.pydlutils.misc.struct_print` test for older Numpy versions.
* Fixed failing yanny file test.
* Improve test infrastructure, including Travis builds.
* Allow comment characters inside quoted strings in yanny files.

0.2.1 (2013-10-06)
------------------

* Added :func:`~pydl.pydlutils.sdss.sdss_sweep_circle`.
* Added first few :mod:`pydl.photoop` functions.
* Clean up some import statements.

0.2.0 (2013-04-22)
------------------

* Using the astropy package-template_ to bring pydl into astropy-compatible form.
* Some but not all tests are re-implemented.

0.1.1 (2013-03-06)
------------------

* Creating a tag representing the state immediately after creation of the
  `git repository`_.

0.1 (2010-11-10)
----------------

* Initial tag (made in svn, not visible in git).  Visible at
  http://www.sdss3.org/svn/repo/pydl/tags/0.1 .

.. _`astropy_helpers`: https://github.com/astropy/astropy-helpers
.. _package-template: https://github.com/astropy/package-template
.. _`git repository`: https://github.com/weaverba137/pydl
.. _`write multiple ndarray to yanny files`: https://github.com/weaverba137/pydl/pull/3
