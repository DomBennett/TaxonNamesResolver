import logging
logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
console.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(console)

from taxon_names_resolver import Resolver
#terms = ['Pt. niger', 'Leis','Harp', 'Not. biguttatus', 'Not. palustris', 'Amara', 'Cal. erratus', 'Pt. diligens', 'Not. reitteri']
terms = ['Harp']
resolver = Resolver(terms = terms, datasource = "NCBI")
resolver.main()