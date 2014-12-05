#! /usr/bin/env python
'''
Example script for using taxon_names_resolver
'''

# SETUP LOGGING
import logging
logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
console.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console)

# PACKAGES
from taxon_names_resolver import Resolver
from taxon_names_resolver import taxTree

# EXAMPLE NAMES
terms = ['Homo sapiens', 'Gorilla gorilla', 'Pongo pongo', 'Macca mulatta',
         'Mus musculus', 'Ailuropoda melanoleuca', 'Ailurus fulgens',
         'Chlorotalpa tytonis', 'Arabidopsis thaliana', 'Bacillus subtilus']

# RESOLVE
resolver = Resolver(terms=terms, datasource="NCBI")
resolver.main()

# CREATE TREE
ranks = resolver.retrieve('classification_path_ranks')
idents = resolver.retrieve('query_name')
lineages = resolver.retrieve('classification_path')
treestring = taxTree(idents, ranks, lineages)

# SAVE
with open('example.tre', 'w') as file:
    file.write(treestring)
