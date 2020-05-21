import sys
import os
try: from setuptools import setup
except ImportError: from distutils.core import setup
import versioneer

URL = 'https://github.com/kadrlica/cartosky'

setup(
    name='cartosky',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url=URL,
    author='Alex Drlica-Wagner',
    author_email='kadrlica@fnal.gov',
    scripts = [],
    install_requires=[
        'matplotlib',
        'numpy',
        'scipy',
        'cartopy >= 0.18.0',
        'ephem',
        'healpy > 1.10.2', # This is required for `lonlat` argument
        'astropy',
        'pandas',
    ],
    packages=['cartosky','cartosky.instrument'],
    package_data={'cartosky':['data/*.txt','data/*.dat']},
    description="Python tools for making skymaps",
    long_description="See %s"%URL,
    platforms='any',
    keywords='python astronomy plotting',
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
    ]
)
