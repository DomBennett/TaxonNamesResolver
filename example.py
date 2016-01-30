#! /usr/bin/env python
# Dom Bennett
# 22/01/2015
# TODO: add tutorial on advanced TaxDict usage
'''
Example script for using taxon_names_resolver
'''
# this is for forward compatibility with python 3
from __future__ import absolute_import
from __future__ import print_function

# SETUP LOGGING (OPTIONAL)
import logging
logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
console.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console)

# PACKAGES
from taxon_names_resolver import Resolver
from taxon_names_resolver import TaxDict
from taxon_names_resolver import taxTree

# EXAMPLE NAMES
terms = ['Homo sapiens', 'Gorilla gorilla', 'Pongo pongo', 'Macca mulatta',
         'Mus musculus', 'Ailuropoda melanoleuca', 'Ailurus fulgens',
         'Chlorotalpa tytonis', 'Arabidopsis thaliana', 'Bacillus subtilus']

# RESOLVE
# pass the terms, the datasource and the logger (optional)
resolver = Resolver(terms=terms, datasource="NCBI", logger=logger)
resolver.main()  # resolve!

# CREATE TAXDICT
# extract the unique names for each term ('idents', query_name is best as it is
#  guaranteed to be unique)
idents = resolver.retrieve('query_name')
# extract the lists of names for all known parental taxonomic groups for each
#  term ('lineages', e.g. Homo, Primate, Mammalia)
lineages = resolver.retrieve('classification_path')
# for Taxonomic IDs instead of names, use:
#  lineages = resolver.retrieve('classification_path_ids')
# extract the lists of corresponding rank names for 'lineages' ('ranks', e.g.
#  species, genus etc.) for each entity
ranks = resolver.retrieve('classification_path_ranks')
# optional extra data slots are also possible, for example a list of 1s and 0s
# it could be anything, just as long as its in the same order
extra = [1, 1, 1, 0, 0, 1, 1, 0, 1, 0]
# create a taxonomy specifying the names and order of 'ranks'. N.B. this is the
#  default and is based on NCBI's taxonomy.
taxonomy = ['subspecies', 'species', 'subgenus', 'genus', 'tribe', 'subfamily',
            'family', 'superfamily', 'parvorder', 'infraorder', 'suborder',
            'order', 'superorder', 'parvclass', 'infraclass', 'subclass',
            'class', 'superclass', 'subphylum', 'phylum', 'kingdom',
            'superkingdom']
# use 'idents', 'ranks', 'lineages' and 'taxonomy' (optional) to construct a
#  TaxDict
taxdict = TaxDict(idents=idents, ranks=ranks, lineages=lineages,
                  taxonomy=taxonomy, extra=extra)

# EXPLORE TAXDICT
# a dictionary for each ident with: 'lineage', 'taxref', 'ident', 'cident' and
#  'rank' (+ 'extra')
# the lineage
taxdict['Homo sapiens']['lineage']  # N.B. not all lineages are named, ''
# the ident is the same format as lineage e.g. it could be an ID
taxdict['Homo sapiens']['ident']
# the 'cident' (Contextual Ident), the highest named taxonomic group unique to
#  this ident among all other idents
taxdict['Arabidopsis thaliana']['cident']  # A. thaliana is the only plant
# the 'taxref', a holder of 'ident' and taxonomic posistion. Requires printing
#  e.g. C. tytonis could only resolved to the genus level (22/01/2015)
print(taxdict['Chlorotalpa tytonis']['taxref'])
# check the taxonomy
print(taxdict.taxonomy)
# check the hierarchy, a dictionary of taxrefs ranked and grouped in the form:
#  {'rank':[([taxref1, taxref2, ....],'lineage1'),
#           ([taxref3, taxref4, ....],'lineage2'), ....]}
print(taxdict.hierarchy)
# we've also added an extra data slot
taxdict['Homo sapiens']['extra']

# CREATE TREE
# use the taxdict to create a Newick string
treestring = taxTree(taxdict)

# SAVE TREE
with open('example.tre', 'w') as file:
    file.write(treestring)
