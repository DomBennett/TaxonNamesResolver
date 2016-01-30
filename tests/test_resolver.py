#! /usr/bin/env python
"""
Tests for Resolver
"""
from __future__ import absolute_import

import unittest
import json
import os
import shutil
import taxon_names_resolver as tnr
import six

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
class dummy_Logger(object):

    def __init__(self):
        pass

    def info(self, msg):
        pass

    def debug(self, msg):
        pass

    def warn(self, msg):
        pass

    def error(self, msg):
        pass


# redefining the unbound search method...
def dummy_search(self, terms, prelim):
    if prelim:
        return first
    else:
        return fourth


# add the datasources to prevent talking to GNR
class Dummy_GnrDataSources(object):
    def __init__(self, logger):
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
        self.logger = dummy_Logger()
        self.resolver1 = tnr.resolver.Resolver(terms=terms, taxon_id='51',
                                               logger=self.logger)
        # second has a store added
        self.resolver2 = tnr.resolver.Resolver(terms=terms, taxon_id='51',
                                               logger=self.logger)
        test_store = tnr.gnr_tools.GnrStore(terms, tax_group='51',
                                            logger=self.logger)
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
        exp_header = 'query_name,classification_path,data_source_title,\
match_type,score,classification_path_ranks,name_string,\
canonical_form,classification_path_ids,prescore,\
data_source_id,taxon_id,gni_uuid'
        exp_results = ['GenusA speciesA,|kingdomA|phylumA|orderA|familyA|genusA|speciesA,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusA speciesA,GenusA speciesA,|51|41|31|21|11|1,3|0|0,1,1,test_number1',
                       'GenusA speciesB,|kingdomA|phylumA|orderA|familyA|genusA|speciesB,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusA speciesB,GenusA speciesB,|51|41|31|21|11|2,3|0|0,1,2,test_number2',
                       'GenusA speciesC,|kingdomA|phylumA|orderA|familyA|genusA|speciesC,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusA speciesC,GenusA speciesC,|51|41|31|21|11|3,3|0|0,1,3,test_number3',
                       'GenusC speciesF,|kingdomA|phylumA|orderA|familyB|genusC|speciesF,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusC speciesF,GenusC speciesF,|51|41|31|22|13|6,3|0|0,1,6,test_number6',
                       'GenusD speciesG,|kingdomA|phylumA|orderB|familyC|genusD|speciesG,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusD speciesG,GenusD speciesG,|51|41|32|23|14|7,3|0|0,1,7,test_number7',
                       'GenusB speciesE,|kingdomA|phylumA|orderA|familyA|genusB|speciesE,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusB speciesE,GenusB speciesE,|51|41|31|21|12|5,3|0|0,1,5,test_number5',
                       'GenusB speciesD,|kingdomA|phylumA|orderA|familyA|genusB|speciesD,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusB speciesD,GenusB speciesD,|51|41|31|21|12|4,3|0|0,1,4,test_number4',
                       'GenusF speciesI,|kingdomA|phylumC|orderD|familyE|genusF|speciesI,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusF speciesI,GenusF speciesI,|51|43|34|25|16|9,3|0|0,1,9,test_number9',
                       'GenusE speciesH,|kingdomA|phylumB|orderC|familyD|genusE|speciesH,\
test_source,1,1.0,|kingdom|phylum|order|family|genus|species,\
GenusE speciesH,GenusE speciesH,|51|42|33|24|15|8,3|0|0,1,8,test_number8']

        self.resolver2.write()
        with open(os.path.join('resolved_names', 'unresolved.txt'), 'r')\
                as f:
            res = f.readline().strip('\n')
        self.assertEqual(res, 'GenusG speciesJ')

        with open(os.path.join('resolved_names', 'search_results.csv'), 'r')\
                as f:
            obs_lines = [line.rstrip() for line in f]
        self.assertEqual(exp_header, obs_lines[0])
        six.assertCountEqual(self, exp_results, obs_lines[1:])

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
        res = list(self.resolver1._store.keys())
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
