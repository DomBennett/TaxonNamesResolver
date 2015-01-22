#! /usr/bin/env python
# D.J. Bennett
# 01/06/2014
# Put on PyPi with http://peterdowns.com/posts/first-time-with-pypi.html
"""
TaxonNamesResolver setup
"""

import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

PACKAGES = find_packages()
PACKAGES = [each for each in PACKAGES if each != 'tests']
PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]

setup(
    name="taxon_names_resolver",
    version="1.0.2",
    author="Dominic John Bennett",
    author_email="dominic.john.bennett@gmail.com",
    description=("Resolve taxonomic names through Global Names Resolver."),
    license="LICENSE.txt",
    keywords="taxonomy biodiversity systematics natural history",
    url="https://github.com/DomBennett/TaxonNamesResolver",
    download_url="https://github.com/DomBennett/TaxonNamesResolver/tarball/1.0.2",
    packages=PACKAGES,
    package_dir=dict(zip(PACKAGES, PACKAGE_DIRS)),
    test_suite='tests',
    scripts=['TaxonNamesResolver.py'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    ],
    install_requires=['setuptools'],
)
