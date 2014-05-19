#! /usr/bin/env python
import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

PACKAGES = find_packages()
PACKAGE_DIRS = [p.replace(".", os.path.sep) for p in PACKAGES]

setup(
    name = "taxon_names_resolver",
    version = "0.0.1",
    author = "Dominic John Bennett",
    author_email = "dominic.john.bennett@gmail.com",
    description = ("Resolve taxonomic names through Global Names Resolver."),
    license = "No license",
    keywords = "taxonomy biodiversity systematics natural history",
    url = "https://github.com/DomBennett/TaxonNamesResolver",
    packages = PACKAGES,
    package_dir = dict(zip (PACKAGES, PACKAGE_DIRS)),
    scripts = ['TaxonNamesResolver.py'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1 - Planning",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 2.7",
    ],
    install_requires=[
          # -*- Extra requirements: -*-
          'setuptools',
      ],
)