#! /usr/bin/env python
# D.J. Bennett
# 16/05/2014
"""
TaxonNamesResolver is a python package for resolving taxonomic
names through Global Names Resolver (Copyright (C) 2012-2013
Marine Biological Laboratory). It was written by Dominic John
Bennett with additional help from Lawrence Hudson.

Copyright (C) 2014  Dominic John Bennett

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

# Import
import os
import sys
import argparse
import logging
import platform
from datetime import datetime
from taxon_names_resolver import Resolver

description = """
TaxonNamesResolver D.J. Bennett (C) 2014

Resolve taxonomic names through the Global Names Resolver
(Copyright (C) 2012-2013 Marine Biological Laboratory) by searching
multiple taxonomic datasources.
"""


def parseArgs():
    """Read arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-names", "-n", help=".txt file of taxonomic names")
    parser.add_argument("-datasource", "-d", help="taxonomic datasource by \
which names will be resolved (default NCBI)")
    parser.add_argument("-taxonid", "-t", help="parent taxonomic ID")
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    return parser


def logSysInfo():
    """Write system info to log file"""
    logging.info('\n' + '#' * 70)
    logging.info(datetime.today().strftime("%A, %d %B %Y %I:%M%p"))
    logging.info('Running on [{0}] [{1}]'.format(platform.node(),
                                                 platform.platform()))
    logging.info('Python [{0}]'.format(sys.version))
    logging.info('#' * 70 + '\n')


def logEndTime():
    """Write end info to log"""
    logging.info('\n' + '#' * 70)
    logging.info('Complete')
    logging.info(datetime.today().strftime("%A, %d %B %Y %I:%M%p"))
    logging.info('#' * 70 + '\n')

if __name__ == '__main__':
    parser = parseArgs()
    args = parser.parse_args()
    if not os.path.isfile(args.names):
        print '[{0}] could not be found!'.format(args.names)
        sys.exit()
    print '\n' + description + '\n'
    if args.datasource:
        datasource = args.datasource
    else:
        datasource = 'NCBI'
    # simple logging, no levels, duplicate to console if verbose
    logfile = 'log.txt'
    logging.basicConfig(filename=logfile, level=logging.INFO,
                        format='%(message)s')
    if args.verbose:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        logging.getLogger('').addHandler(console)
    # log system info
    logSysInfo()
    resolver = Resolver(args.names, datasource, args.taxonid)
    resolver.main()
    resolver.write()
    logEndTime()
    if not args.verbose:
        print '\nComplete\n'
