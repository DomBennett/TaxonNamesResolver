#! /usr/bin/env python
# D.J. Bennett
# 16/05/2014

# Import
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import argparse
import logging
import platform
from datetime import datetime
from taxon_names_resolver import Resolver
from taxon_names_resolver import __version__ as version
from taxon_names_resolver import __doc__ as details

description = """
----------------------------------------------------------------------
TaxonNamesResolver Version {0}, Copyright (C) Bennett 2014
----------------------------------------------------------------------
This program comes with ABSOLUTELY NO WARRANTY. This is free software,
and you are welcome to redistribute it under certain conditions.
For more details, type `TaxonNamesResolver.py --details`.
----------------------------------------------------------------------
""".format(version)


# FUNCTIONS
def parseArgs():
    """Read arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-names", "-n", help=".txt file of taxonomic names")
    parser.add_argument("-datasource", "-d", help="taxonomic datasource by \
which names will be resolved (default NCBI)")
    parser.add_argument("-taxonid", "-t", help="parent taxonomic ID")
    parser.add_argument("--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument('--details', help='display information about the \
program', action='store_true')
    return parser.parse_args()


def logSysInfo():
    """Write system info to log file"""
    logger.info('#' * 70)
    logger.info(datetime.today().strftime("%A, %d %B %Y %I:%M%p"))
    logger.info('Running on [{0}] [{1}]'.format(platform.node(),
                                                platform.platform()))
    logger.info('Python [{0}]'.format(sys.version))
    logger.info('#' * 70 + '\n')


def logEndTime():
    """Write end info to log"""
    logger.info('\n' + '#' * 70)
    logger.info('Complete')
    logger.info(datetime.today().strftime("%A, %d %B %Y %I:%M%p"))
    logger.info('#' * 70 + '\n')

# MAIN
if __name__ == '__main__':
    args = parseArgs()
    if args.details:
        print('\nThis is TaxonNamesResolver, version: [{0}]'.format(version))
        print(details)
        sys.exit()
    if not args.names:
        print('No names file provided!')
        print('Type `TaxonNamesResolver.py -h` for help.')
        sys.exit()
    if not os.path.isfile(args.names):
        print('[{0}] could not be found!'.format(args.names))
        sys.exit()
    print('\n' + description + '\n')
    if args.datasource:
        datasource = args.datasource
    else:
        datasource = 'NCBI'
    # simple logging, no levels, duplicate to console if verbose
    logfile = 'log.txt'
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    loghandler = logging.FileHandler(logfile, 'a')
    loghandler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(loghandler)
    if args.verbose:
        console = logging.StreamHandler()
        console.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(console)
    # log system info
    logSysInfo()
    resolver = Resolver(args.names, datasource, args.taxonid)
    resolver.main()
    resolver.write()
    logEndTime()
    if not args.verbose:
        print('\nComplete\n')
