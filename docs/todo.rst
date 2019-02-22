====
TODO
====

* Increase test coverage.
* :mod:`pydl.pydlutils.mangle` needs more work.

  - Area (solid angle) calculation.
  - Area calculation needs to account for the
    ``use_caps`` attribute.
  - Intersection of caps with other caps.

* Use numpy/scipy Cholesky tools

  - https://trac.sdss.org/browser/repo/sdss/idlutils/trunk/pro/bspline/cholesky_band.pro
  - https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cholesky.html
  - https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.cholesky.html

* Update ``astropy_helpers`` to v2.0.8.
* Check ``groupdim``, ``groupsize`` in :func:`~pydl.pydlutils.math.djs_reject`.
  Make sure integer division works.
* Document :class:`~pydl.pydlutils.bspline.bspline`.
