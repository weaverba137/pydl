====
PyDL
====

|Astropy Status| |License| |Zenodo| |PyPI Status| |Actions Status| |Coveralls Status| |Documentation Status|

.. |Astropy Status| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

.. |License| image:: https://img.shields.io/pypi/l/pydl.svg
    :target: https://pypi.python.org/pypi/pydl
    :alt: License

.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2575873.svg
    :target: https://doi.org/10.5281/zenodo.2575873
    :alt: DOI: 10.5281/zenodo.2575873

.. |PyPI Status| image:: https://img.shields.io/pypi/v/pydl.svg
    :target: https://pypi.python.org/pypi/pydl
    :alt: PyPI Badge

.. |Actions Status| image:: https://github.com/weaverba137/pydl/actions/workflows/ci_workflows.yml/badge.svg
    :target: https://github.com/weaverba137/pydl/actions/workflows/ci_workflows.yml
    :alt: GitHub Actions CI Status

.. |Coveralls Status| image:: https://coveralls.io/repos/weaverba137/pydl/badge.svg?branch=main&service=github
    :target: https://coveralls.io/github/weaverba137/pydl?branch=main
    :alt: Test Coverage Status

.. |Documentation Status| image:: https://readthedocs.org/projects/pydl/badge/?version=latest
    :target: http://pydl.readthedocs.org/en/latest/
    :alt: Documentation Status

Description
-----------

This package consists of Python_ replacements for functions that are part of
the `IDL®`_ built-in library or part of astronomical `IDL®`_ libraries.
The emphasis is on reproducing results of the astronomical library functions.
Only the bare minimum of `IDL®`_ built-in functions are implemented to support this.

There are four astronomical libraries targeted:

* idlutils_ : a general suite of tools heavily used by SDSS_.
* `Goddard utilities`_ : The `IDL®`_ Astronomy User's Libary, maintained by Wayne Landsman and distributed with idlutils_.
* idlspec2d_ : tools for working with SDSS_, BOSS_ and eBOSS_ spectroscopic data.
* photoop_ : tools for working with SDSS_ imaging data.

This package affiliated with the astropy_ project and is registered with PyPI_.

Full Documentation
------------------

Please visit `PyDL on Read the Docs`_

History
-------

This package was initially developed on the SDSS-III_ `svn repository`_.  It was
moved to the new GitHub_ repository on 2013-03-06.  The present location of
the repository is http://github.com/weaverba137/pydl .

License
-------

PyDL is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.

Legal
-----

* IDL is a registered trademark of `NV5 Geospatial`_.

.. _Python: https://www.python.org
.. _`IDL®`: https://www.nv5geospatialsoftware.com/Products/IDL
.. _idlutils: https://www.sdss4.org/dr16/software/idlutils/
.. _SDSS: https://www.sdss.org
.. _`Goddard utilities`: https://asd.gsfc.nasa.gov/archive/idlastro/
.. _idlspec2d: https://svn.sdss.org/public/repo/eboss/idlspec2d/trunk/
.. _BOSS: https://www.sdss4.org/surveys/boss/
.. _eBOSS: https://www.sdss4.org/surveys/eboss/
.. _photoop: https://svn.sdss.org/public/repo/sdss/photoop/trunk/
.. _astropy: http://www.astropy.org
.. _PyPI: https://pypi.python.org/pypi/pydl/
.. _`PyDL on Read the Docs`: https://pydl.readthedocs.io/en/latest/
.. _SDSS-III: http://www.sdss3.org
.. _`svn repository`: https://www.sdss.org/dr16/software/products/
.. _GitHub: https://github.com
.. _`NV5 Geospatial`: https://www.nv5geospatialsoftware.com
