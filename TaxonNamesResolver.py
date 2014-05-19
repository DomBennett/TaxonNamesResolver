#! /usr/bin/env python
## No warranty, no copyright
## Dominic John Bennett
## 16/05/2014
"""
Script for TaxonNamesResolver
"""

# Import
import os,sys
from taxon_names_resolver import Resolver

# Script
print '\n\nHello, this is TaxonNamesResolver! For details please see https://github.com/DomBennett/TaxonNamesResolver\n'
print 'Please give the file of the taxon names to be searched'
input_file = raw_input('File name: ')
if not os.path.isfile(input_file):
	print "[{0}] not a file!".format(input_file)
	sys.exit()
print '\nPlease give the Datasource name from which you\'d ' +\
    'like to resolve, or hit return to use NCBI by default'
datasource = raw_input('Datasource: ')
if datasource == '':
	datasource = 'NCBI'
print '\nPlease give the lowest shared taxonomic group ID ' +\
    'or hit return to skip'
taxon_id = raw_input('Taxon ID: ')
print '\n\nUser input ...'
print 'Input File: ' + input_file + '\nDatasource: ' + datasource +\
    '\nTaxon ID: ' + taxon_id + '\n\nIs this all correct?\n' +\
    'If not, hit Ctrl+C (or Cmd+C for Mac) to exit.'
raw_input('Hit return to continue.')
if taxon_id == '':
	taxon_id = False
resolver = Resolver(input_file, datasource, taxon_id)
print '\n\nProcessing -- in batches ...\n'
resolver.main()
resolver.write()
print '\n\nComplete!\n'