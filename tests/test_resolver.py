#! /usr/bin/env python
"""
Tests for Resolver
"""

import unittest
import json
import os
import shutil
import taxon_names_resolver as tnr

# TEST DATA
# results from the first search
with open(os.path.join(os.path.dirname(__file__), 'data',
          'test_firstsearch.json'), 'r') as file:
    first = json.load(file)
# results from a thrid search on the original datasource with alt name
#  with the supplied name corrected
with open(os.path.join(os.path.dirname(__file__), 'data',
          'test_fourthsearch.json'), 'r') as file:
    fourth = json.load(file)
# results containing multiple records:
# -one with a lower matching score
# -one resolved to a higher taxonomic level
with open(os.path.join(os.path.dirname(__file__), 'data',
          'test_multiple.json'), 'r') as file:
    multiple = json.load(file)

terms = ['GenusA speciesA', 'GenusA speciesB', 'GenusA speciesC',
         'GenusB speciesD', 'GenusB speciesE', 'GenusC speciesF',
         'GenusD speciesG', 'GenusE speciesH', 'GenusF speciesI',
         'GenusG speciesJ']


# STUBS
# redefining the unbound search method...
def dummy_search(self, terms, prelim):
    if prelim:
        return first
    else:
        return fourth


# add the datasources to prevent talking to GNR
class Dummy_GnrDataSources(object):
    def __init__(self):
        pass

    def byName(self, names, invert=False):
        if invert:
            return [1, 2, 3]
        else:
            return [4]


class ResolverTestSuite(unittest.TestCase):

    def setUp(self):
        # patch
        self.true_search = tnr.gnr_tools.GnrResolver.search
        self.True_GnrDataSources = tnr.gnr_tools.GnrDataSources
        tnr.gnr_tools.GnrResolver.search = dummy_search
        tnr.gnr_tools.GnrDataSources = Dummy_GnrDataSources
        # first resolver has no store
        self.resolver1 = tnr.resolver.Resolver(terms=terms, taxon_id='51')
        # second has a store added
        self.resolver2 = tnr.resolver.Resolver(terms=terms, taxon_id='51')
        test_store = tnr.gnr_tools.GnrStore(terms, tax_group='51')
        test_store.add(first)
        self.resolver2._store = test_store

    def tearDown(self):
        # replace patches
        tnr.gnr_tools.GnrResolver.search = self.true_search
        tnr.gnr_tools.GnrDataSources = self.True_GnrDataSources

    def test_resolver_private_count(self):
        # all results will have fewer than 2 records
        res = self.resolver2._count(greater=False, nrecords=2)
        self.assertEqual(len(res), 10)
        # one record has 0 records
        res = self.resolver2._count(greater=True, nrecords=0)
        self.assertEqual(len(res), 9)

    def test_resolver_write(self):
        # write out csv and read in and check
        exp = 'GenusA speciesA,|kingdomA|phylumA|orderA|familyA|genusA|speciesA,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,GenusA \
speciesA,GenusA speciesA,|51|41|31|21|11|1,3|0|0,1,1,test_number1'
        self.resolver2.write()
        with open(os.path.join('resolved_names', 'unresolved.txt'), 'r')\
                as file:
            res = file.readline().strip('\n')
        self.assertEqual(res, 'GenusG speciesJ')
        with open(os.path.join('resolved_names', 'search_results.csv'), 'r')\
                as file:
            _skip = file.readline()
            res = file.readline().strip()
        del _skip
        self.assertEqual(res, exp)
        shutil.rmtree('resolved_names')

    def test_resolver_retrieve(self):
        # retrieve datasource ids, all should be 1
        res = self.resolver2.retrieve('data_source_id')
        res = [each == 1 for each in res]
        self.assertTrue(all(res))

    def test_resolver_main(self):
        # test the searches
        # the results should equal 10
        self.resolver1.main()
        res = self.resolver1._store.keys()
        self.assertEqual(len(res), 10)

    def test_resolver_private_sieve(self):
        # filter the multiple records
        # first replace the multiple records in the store
        self.resolver2._store.replace(multiple)
        multi_records = ['GenusA speciesB', 'GenusF speciesI']
        res = self.resolver2._sieve(multi_records)
        res = [len(each['results']) == 1 for each in res]
        self.assertTrue(all(res))

if __name__ == '__main__':
    unittest.main()
