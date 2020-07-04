import sys
import os
try: from setuptools import setup
except ImportError: from distutils.core import setup
import versioneer

URL = 'https://github.com/kadrlica/cartosky'

with open('requirements.txt') as f:
    install_requires = [req.strip() for req in f.readlines() if req[0] != '#']

setup(
    name='cartosky',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url=URL,
    author='Alex Drlica-Wagner',
    author_email='kadrlica@fnal.gov',
    scripts = [],
    python_requires='>=3.6.0',
    setup_requires=['numpy'],
    install_requires=install_requires,
    packages=['cartosky','cartosky.instrument'],
    package_data={'cartosky':['data/*.txt','data/*.dat']},
    description="Python tools for making skymaps",
    long_description="See %s"%URL,
    platforms='any',
    keywords='python astronomy plotting',
    classifiers = [
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
    ]
)
