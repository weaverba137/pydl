====
PyDL
====

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge


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

.. image:: https://img.shields.io/pypi/v/pydl.svg
    :target: https://pypi.python.org/pypi/pydl
    :alt: PyPI Badge

.. image:: https://img.shields.io/pypi/dm/pydl.svg
    :target: https://pypi.python.org/pypi/pydl
    :alt: PyPI Downloads


Full Documentation
------------------

Please visit `PyDL on Read the Docs`_

.. image:: https://readthedocs.org/projects/pydl/badge/?version=latest
    :target: http://pydl.readthedocs.org/en/latest/
    :alt: Documentation Status


History
-------

This package was initially developed on the SDSS-III_ `svn repository`_.  It was
moved to the new GitHub_ repository on 2013-03-06.  The present location of
the repository is http://github.com/weaverba137/pydl .


Travis Build Status
-------------------

.. image:: https://img.shields.io/travis/weaverba137/pydl.svg
    :target: https://travis-ci.org/weaverba137/pydl
    :alt: Travis Build Status


Test Coverage Status
--------------------

.. image:: https://coveralls.io/repos/weaverba137/pydl/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/weaverba137/pydl?branch=master
    :alt: Test Coverage Status


License
-------
.. image:: https://img.shields.io/pypi/l/pydl.svg
    :target: https://pypi.python.org/pypi/pydl
    :alt: License

PyDL is free software licensed under a 3-clause BSD-style license. For details see
the ``licenses/LICENSE.rst`` file.


Legal
-----

* IDL is a registered trademark of `Exelis Visual Information Solutions`_.

.. _Python: http://python.org
.. _`IDL®`: http://www.exelisvis.com/language/en-us/productsservices/idl.aspx
.. _idlutils: http://www.sdss.org/dr12/software/idlutils/
.. _SDSS: http://www.sdss.org
.. _`Goddard utilities`: http://idlastro.gsfc.nasa.gov/
.. _idlspec2d: https://svn.sdss.org/public/repo/eboss/idlspec2d/trunk/
.. _BOSS: http://www.sdss.org/surveys/boss/
.. _eBOSS: http://www.sdss.org/surveys/eboss/
.. _photoop: https://svn.sdss.org/public/repo/sdss/photoop/trunk/
.. _astropy: http://www.astropy.org
.. _PyPI: https://pypi.python.org/pypi/pydl/
.. _`PyDL on Read the Docs`: http://pydl.readthedocs.org/en/latest/
.. _SDSS-III: http://www.sdss3.org
.. _`svn repository`: http://www.sdss.org/dr12/software/products/
.. _GitHub: http://github.com
.. _`Exelis Visual Information Solutions`: http://www.exelisvis.com
