import logging
logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
console.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console)

from taxon_names_resolver import Resolver
terms = ['Homo sapiens', 'Gorilla gorilla', 'Pongo pongo', 'Macca mulatta',
         'Mus musculus', 'Ailuropoda melanoleuca', 'Ailurus fulgens',
         'Chlorotalpa tytonis', 'Arabidopsis thaliana', 'Bacillus subtilus']
resolver = Resolver(terms=terms, datasource="NCBI")
resolver.main()

from taxon_names_resolver import taxTree
ranks = resolver.retrieve('classification_path_ranks')
idents = resolver.retrieve('query_name')
lineages = resolver.retrieve('classification_path')
taxTree(idents, ranks, lineages, draw=True)
