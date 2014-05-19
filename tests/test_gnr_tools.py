#! /usr/bin/env python
"""
Tests for GNR tools
"""

import unittest,json,os
from taxon_names_resolver import gnr_tools as gt

with open(os.path.join('data','test.json'), 'r') as file:
	test_json = json.load(file)

terms = ['GenusA speciesA', 'GenusA speciesB', 'GenusA speciesC', 'GenusB speciesD',\
'GenusB speciesE','GenusC speciesF','GenusD speciesG','GenusE speciesH','GenusF speciesI',\
'GenusG speciesJ']

class GNRToolsTestSuite(unittest.TestCase):

	def test_datasources(self):
		# create test datasources class
		test_ds = gt.GnrDataSources()
		res_sum = test_ds.summary()
		# assert that we get a list of dictionaries
		self.assertTrue(isinstance(res_sum[0], dict))
		first_ds = res_sum[0]['title']
		# test byNames
		res1 = test_ds.byName(names = first_ds)
		res2 = test_ds.byName(names = first_ds, invert = True)
		self.assertEqual([res1[0],len(res2)], [1, len(res_sum)-1])

	def test_store(self):
		# create test store
		test_store = gt.GnrStore(terms)
		# add test json to store
		test_store.add(test_json)
		res1 = test_store['GenusA speciesA'][0]
		self.assertEqual(res1['taxon_id'], '1')
		# replace first species
		replacement = [{'supplied_name_string':'GenusA speciesA'}]
		test_store.replace(replacement)
		res2 = test_store['GenusA speciesA']
		self.assertEqual(res2, [])

if __name__ == '__main__':
    unittest.main()