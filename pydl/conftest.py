# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
#
# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure.
#
import os

from astropy.version import version as astropy_version
if astropy_version < '3.0':
    #
    # With older versions of Astropy, we actually need to import the pytest
    # plugins themselves in order to make them discoverable by pytest.
    #
    from astropy.tests.pytest_plugins import *
    del pytest_report_header
else:
    # As of Astropy 4.0, the pytest plugins provided by Astropy are
    # now in the pytest-astropy-header package.  This is backward-compatible
    # with Astropy 3.
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS

from astropy.tests.helper import enable_deprecations_as_exceptions

## Uncomment the following line to treat all DeprecationWarnings as
## exceptions. For Astropy v2.0 or later, there are 2 additional keywords,
## as follow (although default should work for most cases).
## To ignore some packages that produce deprecation warnings on import
## (in addition to 'compiler', 'scipy', 'pygments', 'ipykernel', and
## 'setuptools'), add:
##     modules_to_ignore_on_import=['module_1', 'module_2']
## To ignore some specific deprecation warning messages for Python version
## MAJOR.MINOR or later, add:
##     warnings_to_ignore_by_pyver={(MAJOR, MINOR): ['Message to ignore']}
# enable_deprecations_as_exceptions()

def pytest_configure(config):

    config.option.astropy_header = True
    #
    # Customize the following lines to add/remove entries from
    # the list of packages for which version numbers are displayed when running
    # the tests.
    #
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['PyDL'] = 'pydl'
    PYTEST_HEADER_MODULES.pop('Pandas', None)
    PYTEST_HEADER_MODULES.pop('h5py', None)

    from .version import version  #, astropy_helpers_version
    packagename = os.path.basename(os.path.dirname(__file__))
    #
    # Display the version number of the package rather than the version number
    # of Astropy in the top line when running the tests.
    #
    TESTED_VERSIONS[packagename] = version
    # TESTED_VERSIONS['astropy_helpers'] = astropy_helpers_version
