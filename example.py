#! /usr/bin/env python
# Dom Bennett
# 22/01/2015
'''
Example script for using taxon_names_resolver
'''

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
                  taxonomy=taxonomy)

# EXPLORE TAXDICT
# a dictionary for each ident with: its lineage, taxref and cident
# the lineage
taxdict['Homo sapiens']['lineage']  # N.B. not all lineages are named, ''
# the 'cident' (Contextual Ident), the highest named taxonomic group unique to
#  this ident among all other idents
taxdict['Arabidopsis thaliana']['cident']  # A. thaliana is the only plant
# the 'taxref', a holder of 'ident' and taxonomic posistion. Requires printing
#  e.g. C. tytonis could only resolved to the genus level (22/01/2015)
print taxdict['Chlorotalpa tytonis']['taxref']
# check the taxonomy
print taxdict.taxonomy
# check the hierarchy, a dictionary of taxrefs ranked and grouped in the form:
#  {'rank':[([taxref1, taxref2, ....],'lineage1'),
#           ([taxref3, taxref4, ....],'lineage2'), ....]}
print taxdict.hierarchy

# CREATE TREE
# use the taxdict to create a Newick string
treestring = taxTree(taxdict)

# SAVE TREE
with open('example.tre', 'w') as file:
    file.write(treestring)
