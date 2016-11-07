#! /usr/bin/env python
# D.J. Bennett
# 01/06/2014
# Put on PyPi with http://peterdowns.com/posts/first-time-with-pypi.html
"""
TaxonNamesResolver setup
"""
from __future__ import absolute_import

# PACKAGES
import os
import taxon_names_resolver
from setuptools import setup, find_packages
from six.moves import zip

# FIND
PACKAGES = find_packages()
PACKAGES = [each for each in PACKAGES if each != 'tests']
PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]
url = 'https://github.com/DomBennett/TaxonNamesResolver/tarball/' + \
    str(taxon_names_resolver.__version__)

# SETUP
setup(
    name="taxon_names_resolver",
    version=taxon_names_resolver.__version__,
    author="Dominic John Bennett",
    author_email="dominic.john.bennett@gmail.com",
    description=("Resolve taxonomic names through Global Names Resolver."),
    license="LICENSE.txt",
    keywords="taxonomy biodiversity systematics natural history",
    url="https://github.com/DomBennett/TaxonNamesResolver",
    download_url=url,
    packages=PACKAGES,
    package_dir=dict(list(zip(PACKAGES, PACKAGE_DIRS))),
    test_suite='tests',
    scripts=['TaxonNamesResolver.py'],
    long_description=taxon_names_resolver.__doc__,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    ],
    install_requires=['setuptools', 'six'],
)
