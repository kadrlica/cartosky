# cartosky

[![Build](https://img.shields.io/travis/kadrlica/cartosky.svg)](https://travis-ci.org/kadrlica/cartosky)
[![PyPI](https://img.shields.io/pypi/v/cartosky.svg)](https://pypi.python.org/pypi/cartosky)
[![Release](https://img.shields.io/github/release/kadrlica/cartosky.svg)](../../releases)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](../../)

The `cartosky` package provides a astronomically oriented interface to ploting sky maps based on [`matplotlib.basemap`](http://matplotlib.org/basemap/). This package addresses several issues present in the [`healpy`](https://healpy.readthedocs.io/en/latest/) plotting routines:
1. `healpy` supports a limited set of sky projections (`cartview`, `mollview`, and `gnomview`)
2. `healpy` converts sparse healpix maps to full maps to plot; this is memory intensive for large `nside`

In addition, `cartosky` provides some convenience functionality for large optical surveys.

## Installation

The best way to install cartosky is if you have [anaconda](https://anaconda.org/) installed. If you have trouble, check out the [.travis.yml](.travis.yml) file. The procedure below will create a conda environment and pip install cartosky:
```
conda create -n cartosky numpy scipy pandas matplotlib cartopy astropy ephem healpy nose -c conda-forge
source activate cartosky
pip install cartosky
```
If you want the bleeding edge of cartosky, you can follow the directions above to create the conda environment, but then install by cloning directly from github:
```
git clone https://github.com/kadrlica/cartosky.git
cd cartosky
python setup.py install
```

## Tutorial

If you want to see what you can do with `cartosky`, check out the [tutorial](tutorial/).
