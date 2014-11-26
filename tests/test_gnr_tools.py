#! /usr/bin/env python
"""
Tests for GNR tools
"""

import unittest
import json
import os
from taxon_names_resolver import gnr_tools as gt

# TEST DATA
# results from the first search
with open(os.path.join(os.path.dirname(__file__), 'data',
          'test_firstsearch.json'), 'r') as file:
    first = json.load(file)
# results from the secondary search in other datasources for alt names
with open(os.path.join(os.path.dirname(__file__), 'data',
          'test_secondsearch.json'), 'r') as file:
    second = json.load(file)
# results from a thrid search on the original datasource with alt name
with open(os.path.join(os.path.dirname(__file__), 'data',
          'test_thirdsearch.json'), 'r') as file:
    third = json.load(file)

terms = ['GenusA speciesA', 'GenusA speciesB', 'GenusA speciesC',
         'GenusB speciesD', 'GenusB speciesE', 'GenusC speciesF',
         'GenusD speciesG', 'GenusE speciesH', 'GenusF speciesI',
         'GenusG speciesJ']

# raw output from GNR
real_json = {u'status': u'success', u'parameters': {u'header_only': False, u'best_match_only':\
False, u'data_sources': [4], u'with_context': False, u'preferred_data_sources': [],\
u'resolve_once': False}, u'data_sources': [{u'id': 4, u'title': u'NCBI'}], u'url':\
u'http://resolver.globalnames.org/name_resolvers/14uh5zifltli.json', u'message': u'Success',\
u'data': [{u'supplied_name_string': u'Homo sapiens', u'results': [{u'classification_path':\
u'|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Coelomata|Deuterostomia|Chordata|Craniata\
|Vertebrata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Tetrapoda|Amniota|Mammalia|\
Theria|Eutheria|Euarchontoglires|Primates|Haplorrhini|Simiiformes|Catarrhini|Hominoidea|Hominidae\
|Homininae|Homo|Homo sapiens', u'data_source_title': u'NCBI', u'match_type': 1, u'score': 0.988,\
u'classification_path_ranks': u'|superkingdom||kingdom|||||phylum|subphylum||superclass||||||class\
|||superorder|order|suborder|infraorder|parvorder|superfamily|family|subfamily|genus|species',\
u'name_string': u'Homo sapiens', u'canonical_form': u'Homo sapiens', u'classification_path_ids':\
u'131567|2759|33154|33208|6072|33213|33316|33511|7711|89593|7742|7776|117570|117571|8287|32523|\
32524|40674|32525|9347|314146|9443|376913|314293|9526|314295|9604|207598|9605|9606', u'prescore':\
u'3|0|0', u'data_source_id': 4, u'taxon_id': u'9606', u'gni_uuid': u'16f235a0-e4a3-529c-9b83-bd15\
fe722110'}]}, {u'supplied_name_string': u'Gallus gallus', u'results': [{u'classification_path': u'\
|Eukaryota|Opisthokonta|Metazoa|Eumetazoa|Bilateria|Coelomata|Deuterostomia|Chordata|Craniata|Verte\
brata|Gnathostomata|Teleostomi|Euteleostomi|Sarcopterygii|Tetrapoda|Amniota|Sauropsida|Sauria|Archo\
sauria|Dinosauria|Saurischia|Theropoda|Coelurosauria|Aves|Neognathae|Galliformes|Phasianidae|Phasia\
ninae|Gallus|Gallus gallus', u'data_source_title': u'NCBI', u'match_type': 1, u'score': 0.988, u'cl\
assification_path_ranks': u'|superkingdom||kingdom|||||phylum|subphylum||superclass|||||||||||||cl\
ass|superorder|order|family|subfamily|genus|species', u'name_string': u'Gallus gallus', u'canonical\
_form': u'Gallus gallus', u'classification_path_ids': u'131567|2759|33154|33208|6072|33213|33316|335\
11|7711|89593|7742|7776|117570|117571|8287|32523|32524|8457|32561|8492|436486|436489|436491|436492|8\
782|8825|8976|9005|9072|9030|9031', u'prescore': u'3|0|0', u'data_source_id': 4, u'taxon_id': u'9031\
', u'gni_uuid': u'6b84dfff-e38b-5560-be91-35e423f9ff2a'}]}], u'id': u'14uh5zifltli'}


# STUBS
class Dummy_GnrResolver(gt.GnrResolver):
    pass


def dummy_query(self, terms, data_source_ids):
    # return  real json dictionary
    return real_json

Dummy_GnrResolver._query = dummy_query


class GNRToolsTestSuite(unittest.TestCase):
    # no tests for search and write

    def setUp(self):
        self.resolver = gt.GnrResolver()
        self.dummy_resolver = Dummy_GnrResolver()

    def test_datasources(self):
        # create test datasources class
        test_ds = gt.GnrDataSources()
        res_sum = test_ds.summary()
        # assert that we get a list of dictionaries
        self.assertTrue(isinstance(res_sum[0], dict))
        first_ds = res_sum[0]['title']
        # test byNames
        res1 = test_ds.byName(names=first_ds)
        res2 = test_ds.byName(names=first_ds, invert=True)
        self.assertEqual([res1[0], len(res2)], [1, len(res_sum)-1])

    def test_store(self):
        # create test store
        test_store = gt.GnrStore(terms)
        # add test json to store
        test_store.add(first)
        res1 = test_store['GenusA speciesA'][0]
        self.assertEqual(res1['taxon_id'], '1')
        # replace first species
        replacement = [{'supplied_name_string': 'GenusA speciesA'}]
        test_store.replace(replacement)
        res2 = test_store['GenusA speciesA']
        self.assertEqual(res2, [])

    def test_resolver_private_query(self):
        # talk to gnr with small query
        # Humans are 9606 in NCBI
        res = self.resolver._query(['Homo sapiens'], [4])
        txid = res['data'][0]['results'][0]['taxon_id']
        self.assertEqual(txid, '9606')

    def test_resolver_private_resolve(self):
        # give resolve a list, will return dummy_query results
        res = self.dummy_resolver._resolve([1], [1])
        res = res[0]['supplied_name_string']
        self.assertEqual(res, 'Homo sapiens')

    def test_resolver_private_parsename(self):
        # this function finds names that are different
        #  from the supplied name
        # test_secondsource_json has a different name to supplied
        res = self.resolver._parseNames(second)
        self.assertEqual(res, [(u'GenusG speciesJ', u'GenusG speciesJbis')])

    def test_resolver_private_replacesupstrname(self):
        # this function replaces the new name from a secondary datasource
        #  in the supplied name slot to the original name provided by
        #  the user
        alt_terms = [(u'GenusG speciesJ', u'GenusG speciesJbis')]
        res = self.resolver._replaceSupStrNames(third, alt_terms)
        self.assertEqual(res[0]['supplied_name_string'], 'GenusG speciesJ')

if __name__ == '__main__':
    unittest.main()
