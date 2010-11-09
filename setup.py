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
)

