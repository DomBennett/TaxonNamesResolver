#!/usr/bin/python
## No warranty, no copyright
## Dominic John Bennett
## 07/06/2013
## TODO: use a minimal score (avoid calling hybrids 'species')
## TODO: create a stats-out function
import contextlib
import json
import os
import csv
import urllib, urllib2
import re
import copy

class GnrDataSources(object):
	"""GNR data sources class: extract IDs for specified data sources."""
	def __init__(self):
		url = 'http://resolver.globalnames.org/data_sources.json'
		with contextlib.closing(urllib2.urlopen(url)) as f:
			res = f.read()
		self.available = json.loads(res)

	def summary(self):
		return [dict(id=ds['id'], title=ds['title']) for ds in self.available]

	def byName(self, names, invert = False):
		if invert:
			return [ds['id'] for ds in self.available if not ds['title'] in names]
		else:
			return [ds['id'] for ds in self.available if ds['title'] in names]
			
class GnrResolver(object):
	"""GNR resolver class: search the GNR"""
	def __init__(self, datasource = 'NCBI'):
		ds = GnrDataSources()
		self.write_counter = 1
		self.Id = ds.byName(datasource)
		self.otherIds = ds.byName(datasource, invert = True)

	def search(self, terms, prelim = True):
		"""Search terms against GNR. If prelim = False, search other datasources for alternative names (i.e. synonyms) with which to search main datasource. Return JSON object."""
		if prelim: # preliminary search
			res = self._resolve(terms, self.Id)
			self._write(res)
			return res
		else: # search other DSs for alt names, search DS with these
			res = self._resolve(terms, self.otherIds)
			self._write(res)
			alt_terms = self._parseNames(res)
			if len(alt_terms) == 0:
				return False
			else:
				terms = [each[1] for each in alt_terms] # unzip
				res = self._resolve(terms, self.Id)
				self._write(res)
				alt_res = self._replaceSupStrNames(res, alt_terms)
				return alt_res
	
	def _parseNames(self, jobj):
		# return a list of tuples (term, name) from second search
		# TODO(07/06/2013): record DSs used 
		alt_terms = []
		for record in jobj:
			if len(record) < 2:
				pass
			else:
				term = record['supplied_name_string']
				results = record['results']
				for result in results:
					r_name = result['canonical_form']
					if r_name == term:
						pass
					else:
						alt_terms.append((term, r_name))
		alt_terms = list(set(alt_terms))
		return alt_terms

	def _replaceSupStrNames(self, jobj, alt_terms):
		# replace sup name in jobj with original terms
		for record in jobj:
			sup_name = record['supplied_name_string']
			term = [i for i, each in enumerate(alt_terms) if each[1] == sup_name]
			# avoid the possibility of having the same term with >1 r_names
			term = alt_terms.pop(term[0])[0]
			record['supplied_name_string'] = term
		return jobj
		

	def _resolve(self, terms, ds_id):
		# Query server in chunks
		chunk_size = 100
		res = []
		lower = 0
		while lower < len(terms):
			upper = min(len(terms), lower + chunk_size)
			print 'Querying [{0}] to [{1}] of [{2}]'.format(lower, upper, len(terms))
			res.append(self._query(terms[lower:upper], ds_id))
			lower = upper
		res = [record for search in res for record in search['data']]
		return(res)		

	def _query(self, terms, data_source_ids):
		ds_ids = [str(id) for id in data_source_ids]
		terms = [urllib.quote(unicode(t).encode('utf8')) for t in terms]
		url = ('http://resolver.globalnames.org/name_resolvers.json?' + 
		'data_source_ids=' + '|'.join(ds_ids) + '&' + 
		'resolve_once=false&' + 
		'names=' + '|'.join(terms))
		with contextlib.closing(urllib2.urlopen(url)) as f:
			return json.loads(f.read())
	
	def _write(self, jobj):
		directory = os.path.join(os.getcwd(), 'resolved_names')
		filename = "{0}_raw_results.json".format(self.write_counter)
		jobj_file = os.path.join(directory, filename)
		with open(jobj_file, 'w') as outfile:
			json.dump(jobj, outfile)
		self.write_counter += 1
		

class GnrStore(dict):
	"""GNR store class: acts like a dictionary for GNR JSON format"""
	def __init__(self, terms):
		for term in terms:
			self[term] = []
		
	def add(self, jobj):
		if not isinstance(jobj, bool):
			for record in jobj:
				term = record['supplied_name_string']
				try:
					if len(record) > 1:
						self[term].extend(record['results'])
					else:
						self[term] = []
				except KeyError:
					print 'JSON object contains terms not in GnrStore'
	
	def replace(self, jobj):
		for record in jobj:
			term = record['supplied_name_string']
			try:
				if len(record) > 1:
					self[term] = record['results']
				else:
					self[term] = []
			except KeyError:
				print 'JSON object contains terms not in GnrStore'
		
class TaxonNamesResolver(object):
	"""Taxon Names Resovler class : Automatically resolves taxon names through GNR. All output written in 'resolved_names' folder. See https://github.com/DomBennett/TaxonNamesResolver for details."""
	def __init__(self, input_file = False, datasource = False, \
	 taxon_id = False):
		# organising dirs
		self.directory = os.getcwd()
		self.outdir = os.path.join(self.directory, 'resolved_names')
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)
		input_file = os.path.join(self.directory, input_file)
		# reading in terms
		terms = []
		with open(input_file) as names:
			for name in names:
				terms.append(name.strip())
		terms = [term for term in terms if not term == '']
		print '\nFound [{0}] taxon names to search in input file... '.format(len(terms))
		terms = list(set(terms))
		print '... of which [{0}] are unique.'.format(len(terms))
		# init dep classes
		self.terms = terms
		self._res = GnrResolver(datasource)
		self.primary_datasource = datasource
		self._store = GnrStore(terms)
		self.taxon_id = taxon_id
		self.key_terms = ['query_name', 'classification_path', 'data_source_title', \
					  'match_type', 'score', 'classification_path_ranks',\
					  'name_string', 'canonical_form',\
					  'classification_path_ids', 'prescore','data_source_id',\
					  'taxon_id', 'gni_uuid'] # http://resolver.globalnames.org/api
		# self.tnr_obj = [] # this will hold all output
		
	def main(self):
		"""Search and sieve query names."""
		primary_bool = True
		no_records = True
		nsearch = 1
		search_terms = self.terms
		while no_records:
			if primary_bool:
				print 'Searching [{0}] ...'.format(self.primary_datasource)
			else:
				print 'Searching other datasources ...'
			res = self._res.search(search_terms, prelim = primary_bool)
			if nsearch > 2 and res:
				for each_res, original_name in zip(res, original_names):
					each_res['supplied_name_string'] = original_name
			self._store.add(res)
			no_records = self._count(nrecords = 1) # Check for returns without records
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
			nsearch += 1
			search_terms = no_records
		# Check for multiple records
		multi_records = self._count(greater = True, nrecords = 1)
		if multi_records:
			print 'Choosing best records to return ...'
			res = self._sieve(multi_records, self.taxon_id)
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
	
	def _sieve(self, multiple_records, tax_group = False):
		"""Return json object without multiple returns per resolved name. Names with multiple records are reduced by finding the name in the clade of interest, have the highest score, have the lowest taxonomic rank and/or are the first item returned."""
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
				# in correct taxonomic group?
				if tax_group:
					key_str = 'classification_path_ids'
					bool_tg = [0] * len(results)
					for i,result in enumerate(results):
						class_ids = result[key_str].split('|')
						class_ids = [int(tid) for tid in class_ids]
						if tax_group in class_ids:
							bool_tg[i] = 1
					results = boolResults(results, bool_tg)
				# choose result with best score
				scores = [result['score'] for result in results]
				bool_score = [1 if score == max(scores) else 0 for score in scores]
				results = boolResults(results, bool_score)
				# choose result resolved to lowest taxonomic rank
				res_ranks = [result['classification_path_ranks'].split('|')[-1] for result in results]
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
		"""Return data for key term specified for each resolved name as a list. Possible terms (02/12/2013): 'query_name', 'classification_path', 'data_source_title', 'match_type', 'score', 'classification_path_ranks', 'name_string', 'canonical_form', 'classification_path_ids', 'prescore', 'data_source_id', 'taxon_id', 'gni_uuid'"""
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
			retrieved = [[r2 for r2 in r1.split('|')] for r1 in retrieved]
		return retrieved

if __name__ == "__main__":
	print '\n\nHello, this is TaxonNamesResolver! For details please see https://github.com/DomBennett/TaxonNamesResolver\n'
	print 'Please give the file of the taxon names to be searched'
	input_file = raw_input('File name: ')
	#input_file = "0_data/mammal_taxnames.txt"
	print '\nPlease give the Datasource name from which you\'d ' +\
	    'like to resolve, or hit return to use NCBI by default'
	datasource = raw_input('Datasource: ')
	#datasource = ""
	if datasource == '':
		datasource = 'NCBI'
	print '\nPlease give the lowest shared taxonomic group ID ' +\
	    'or hit return to skip'
	taxon_id = raw_input('Taxon ID: ')
	#taxon_id = "40674"
	print '\n\nUser input ...'
	print 'Input File: ' + input_file + '\nDatasource: ' + datasource +\
	    '\nTaxon ID: ' + taxon_id + '\n\nIs this all correct?\n' +\
	    'If not, hit Ctrl+C (or Cmd+C for Mac) to exit.'
	raw_input('Hit return to continue.')
	if taxon_id == '':
		taxon_id = False
	resolver = TaxonNamesResolver(input_file, datasource, taxon_id)
	print '\n\nProcessing (N.B. queries are in batches) ...\n'
	resolver.main()
	resolver.write()
	print '\n\nComplete!\n'
