#! /usr/bin/env python
## Dominic John Bennett
## 16/05/2014
## TODO: use a minimal score (avoid calling hybrids 'species')
## TODO: create a stats-out function
"""
Resolver class for parsing GNR records.
"""

import json,os,csv,re,copy,logging
from gnr_tools import GnrStore
from gnr_tools import GnrResolver
import urllib

class EncodingError(Exception):
	pass

class Resolver(object):
	"""Taxon Names Resovler class : Automatically resolves taxon names \
through GNR. All output written in 'resolved_names' folder.
See https://github.com/DomBennett/TaxonNamesResolver for details."""
	def __init__(self, input_file = None, datasource = 'NCBI', \
	 taxon_id = None, terms = None):
		# organising dirs
		self.directory = os.getcwd()
		self.outdir = os.path.join(self.directory, 'resolved_names')
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)
		if input_file:
			input_file = os.path.join(self.directory, input_file)
			# reading in terms
			terms = []
			with open(input_file) as names:
				for name in names:
					terms.append(name.strip())
			terms = [term for term in terms if not term == '']
		else:
			if not terms:
				logging.info("No terms provided")
		terms = list(set(terms))
		logging.info('\nFound [{0}] taxon names to search in input file... '.\
		format(len(terms)))
		logging.info('... of which [{0}] are unique.'.format(len(terms)))
		# init dep classes
		self._check(terms)
		self.terms = terms
		self._res = GnrResolver(datasource)
		self.primary_datasource = datasource
		self._store = GnrStore(terms, tax_group = taxon_id)
		self.key_terms = ['query_name', 'classification_path', 'data_source_title', \
					  'match_type', 'score', 'classification_path_ranks',\
					  'name_string', 'canonical_form',\
					  'classification_path_ids', 'prescore','data_source_id',\
					  'taxon_id', 'gni_uuid'] # http://resolver.globalnames.org/api
		# self.tnr_obj = [] # this will hold all output

	def _check(self, terms):
		"""Check terms do not contain unknown characters"""
		for t in terms:
			try:
				_ = urllib.quote(unicode(t).encode('utf8'))
			except:
				logging.error('Unknown character in [{0}]!'.format(t))
				logging.error('.... remove character and try again.')
				raise EncodingError
		
	def main(self):
		"""Search and sieve query names."""
		primary_bool = True
		no_records = True
		nsearch = 1
		search_terms = self.terms
		original_names = []
		while True:
			if primary_bool:
				logging.info('Searching [{0}] ...'.format(self.primary_datasource))
			else:
				logging.info('Searching other datasources ...')
			res = self._res.search(search_terms, prelim = primary_bool)
			if nsearch > 2 and res:
				# if second search failed, look up alternative names
				for each_res, original_name in zip(res, original_names):
					each_res['supplied_name_string'] = original_name
			self._store.add(res)
			no_records = self._count(nrecords = 1) # Check for returns without records
			if no_records:
				if nsearch == 1:
					primary_bool = False
				elif nsearch == 2:
					original_names = no_records
					no_records = [e.split()[0] for e in no_records]  # genus names
					primary_bool = True
				elif nsearch == 3:
					original_names = no_records
					no_records = [e.split()[0] for e in no_records]
					primary_bool = False
				else:
					break
			else:
				break
			nsearch += 1
			search_terms = no_records
		# Check for multiple records
		multi_records = self._count(greater = True, nrecords = 1)
		if multi_records:
			logging.info('Choosing best records to return ...')
			res = self._sieve(multi_records)
			self._store.replace(res)
		
	#def extract(self, what): # depends on tnr
	#	lkeys = ['qnames', 'rnames', 'taxonids', 'ranks']
	#	i = [i for i, each in enumerate(lkeys) if what is each][0]
	#	res = [each[i] for each in self.tnr_obj]
	#	return res
			
	def _readInJson(self):
		jobj_file = os.path.join(self.outdir, 'prelim_search_results.json')
		with open(jobj_file, 'r') as infile:
			jobj = json.load(infile)
		return jobj

	def _count(self, greater = False, nrecords = 0):
		# return a list of all keys that have records 
		#  greater or less than nrecords
		GnrStore = self._store
		assessed = []
		lens = [len(GnrStore[key]) for key in GnrStore.keys()]
		if greater:
			len_bools = [each > nrecords for each in lens]
		else:
			len_bools = [each < nrecords for each in lens]
		for i, key in enumerate(GnrStore.keys()):
			if len_bools[i]:
				assessed.append(key)
		if len(assessed) == 0:
			return False
		else:
			return assessed
	
	def _sieve(self, multiple_records):
		"""Return json object without multiple returns per resolved name.\
Names with multiple records are reduced by finding the name in the clade of\ 
interest, have the highest score, have the lowest taxonomic rank and/or are\
the first item returned."""
		GnrStore = self._store
		def writeAsJson(term, results):
			record = {'supplied_name_string': term}
			if len(results) > 0:
				record['results'] = results
			return record
		def boolResults(results, bool_li, rand = False):
			if rand:
				results = [results[0]] # choose first record (most likely best?)
			elif sum(bool_li) == 1:
				results = [results[bool_li.index(1)]]
			elif sum(bool_li) == 0:
				return [] # return 'no_record'
			else:
				results = [result for i, result in enumerate(results) if bool_li[i]]
			return results
		sieved = []
		ranks = ['species', 'genus', 'family', 'order', 'superorder', 'class',\
		'superclass', 'subphylum', 'phylum', 'kingdom', 'superkingdom']
		for term in multiple_records:
			results = GnrStore[term]
			while len(results) > 1:
				# choose result with best score
				scores = [result['score'] for result in results]
				bool_score = [1 if score == max(scores) else 0 for score in scores]
				results = boolResults(results, bool_score)
				# choose result resolved to lowest taxonomic rank
				res_ranks = [result['classification_path_ranks'].split('|')[-1] for\
				result in results]
				for j, rank in enumerate(ranks):
					bool_rank = [1 if res_rank == rank else 0 for res_rank in res_ranks]
					if sum(bool_rank) > 0:
						break
				results = boolResults(results, bool_rank)
				results = boolResults(results, bool_rank, rand = True)
			record = writeAsJson(term, results)
			sieved.append(record)
		return sieved
		
	#def _parseResults(self):
	# Deprecated (02/12/2013): unnecessary, shouldn't determine the output for the user
	#	def buildTnrObj(GnrStore, tnr_obj, genera = False):
	#		keys = GnrStore.keys()
	#		for key in keys:
	#			results = GnrStore[key]
	#			if len(results) == 0:
	#				pass
	#			else:
	#				results = results[0]
	#				if genera:
	#					ranks = results['classification_path_ranks']
	#					rank = ranks.split('|')[-1]
	#					if rank not in genera:
	#						rname = results['canonical_form']
	#						taxonid = results['taxon_id']
	#						tnr_obj.append([key, rname, taxonid, rank])
	#						genera.append(rname)
	#						del GnrStore[key]
	#				else:
	#					ranks = results['classification_path_ranks']
	#					rank = ranks.split('|')[-1]
	#					if rank == 'species':
	#						rname = results['canonical_form']
	#						taxonid = results['taxon_id']
	#						tnr_obj.append([key, rname, taxonid, rank])
	#						del GnrStore[key]
	#		return tnr_obj
	#	GnrStore = self._store.copy()
	#	tnr_obj = []
	#	tnr_obj = buildTnrObj(GnrStore, tnr_obj)
	#	genera = [each[1] for each in tnr_obj]
	#	genera = [re.split('\s', each)[0] for each in genera]
	#	tnr_obj = buildTnrObj(GnrStore, tnr_obj, genera)
	#	self.tnr_obj.extend(tnr_obj)
		
	def write(self):
		"""Write csv file of resolved names and txt file of unresolved names."""
		csv_file = os.path.join(self.outdir, 'search_results.csv')
		txt_file = os.path.join(self.outdir, 'unresolved.txt')
		headers = self.key_terms
		unresolved = []
		with open(csv_file, 'wb') as file:
			writer = csv.writer(file)
			writer.writerow(headers)
			for key in self._store.keys():
				results = self._store[key]
				if len(results) == 0:
					unresolved.append(key)
				else:
					row = [key]
					for key_term in headers[1:]:
						element = results[0][key_term]
						# GNR returns UTF-8, csv requires ascii
						if 'encode' in dir(element):
							element = element.encode('ascii')
						row.append(element)
					writer.writerow(row)
		if len(unresolved) > 0:
			with open(txt_file, 'w') as file:
				for name in unresolved:
					file.write("{0}\n".format(name))

	def retrieve(self, key_term):
		"""Return data for key term specified for each resolved name as a list.\
Possible terms (02/12/2013): 'query_name', 'classification_path', 'data_source_title',\
'match_type', 'score', 'classification_path_ranks', 'name_string', 'canonical_form',\
'classification_path_ids', 'prescore', 'data_source_id', 'taxon_id', 'gni_uuid'"""
		if key_term not in self.key_terms:
			raise IndexError('Term given is invalid! Check doc string for valid terms.')
		store = copy.deepcopy(self._store)
		retrieved = []
		for key in store.keys():
			record = store[key]
			if len(record) > 0:
				if key_term == 'query_name':
					retrieved.append(key)
				else:
					retrieved.append(record[0][key_term])
		if re.search('path', key_term):
			retrieved = [[r2 for r2 in r1.split('|')[1:]] for r1 in retrieved]
		return retrieved