from distutils.core import setup

setup(
    name='pydl',
    version='0.1',
    author='Benjamin Alan Weaver',
    author_email='benjamin.weaver@nyu.edu',
    packages=['pydl', 'pydlutils', 'pydlspec2d', ],
    url='http://pypi.python.org/pypi/pydl/',
    license='LICENSE.txt',
    description='Replacement for IDL',
    long_description=open('README.txt').read(),
    classifiers = [ 'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
)

