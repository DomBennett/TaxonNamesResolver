#!/usr/bin/python
## No warranty, no copyright
## Dominic John Bennett
## 03/12/2013
## TODO: adapt genTaxTree to use taxonomic information outside of named ranks

import collections, re, copy
from Bio import Phylo
from cStringIO import StringIO

class TestTaxonNamesResolver(object):
	"""Test Taxon Names Resovler class : for doctesting taxon names resolver functions."""
	def __init__(self, scenario = 1):
		if scenario == 1: # simple scenario
			query_names = ["term1", "term2", "term3", "term4", "term5"]
			name_strings = ["sp1", "sp2", "sp3", "sp4", "sp5"]
			taxon_ids = [1,2,3,4,5]
			classification_paths = [["claA", "ordA", "famA", "genA", "spA"],\
				["claA", "ordA", "famA", "genA", "spB"],\
				["claA", "ordA", "famA", "genB", "spC"],\
				["claA", "ordA", "famB", "genC", "spD"],\
				["claA", "ordB", "famC", "genD", "spE"]]
			classification_path_ranks = [["class", "order", "family", "genus", "species"],\
				["class", "order", "family", "genus", "species"],\
				["class", "order", "family", "genus", "species"],\
				["class", "order", "family", "genus", "species"],\
				["class", "order", "family", "genus", "species"]]
			classification_path_ids = [["11", "21", "31", "41", "1"],\
				["11", "21", "31", "41", "2"],\
				["11", "21", "31", "42", "3"],\
				["11", "21", "32", "43", "4"],\
				["11", "22", "33", "44", "5"]]
		if scenario == 2: # low resolution for one term
			query_names = ["term1", "term2", "term3", "term4", "term5"]
			name_strings = ["gen1", "sp2", "sp3", "sp4", "sp5"]
			taxon_ids = [41,2,3,4,5]
			classification_paths = [["claA", "ordA", "famA", "genA"],\
				["claA", "ordA", "famA", "genA", "spB"],\
				["claA", "ordA", "famA", "genB", "spC"],\
				["claA", "ordA", "famB", "genC", "spD"],\
				["claA", "ordB", "famC", "genD", "spE"]]
			classification_path_ranks = [["class", "order", "family", "genus"],\
				["class", "order", "family", "genus", "species"],\
				["class", "order", "family", "genus", "species"],\
				["class", "order", "family", "genus", "species"],\
				["class", "order", "family", "genus", "species"]]
			classification_path_ids = [["11", "21", "31", "41"],\
				["11", "21", "31", "41", "2"],\
				["11", "21", "31", "42", "3"],\
				["11", "21", "32", "43", "4"],\
				["11", "22", "33", "44", "5"]]
		store_obj = {}
		for i in range(len(query_names)):
			store_obj[query_names[i]] = {'name_string':name_strings[i], 'taxon_id':taxon_ids[i],\
				"classification_path":classification_paths[i], 'classification_path_ranks':classification_path_ranks[i],\
				'classification_path_ids':classification_path_ids[i]}
		self._store = store_obj
		self.key_terms = ['query_name', 'classification_path','classification_path_ranks',\
					  'name_string','classification_path_ids','taxon_id']

	def retrieve(self, key_term):
		if key_term not in self.key_terms:
			raise IndexError('Term given is invalid! Check doc string for valid terms.')
		store = copy.deepcopy(self._store)
		retrieved = []
		for key in store.keys():
			if key_term == 'query_name':
				retrieved.append(key)
			else:
				record = store[key]
				retrieved.append(record[key_term])
		return retrieved


def genTaxTree(resolver, by = 'qnames', draw = False):
	"""Generate Newick tree from TaxonNamesResolver class.
		
		Arguments:
		 resolver = TaxonNamesResolver class
		 by = Tip labels, either 'qnames', 'taxids' or 'name_string'
		 draw = Draw ascii tree (logical)

		 Return:
		  (Newick Tree Object, [shared lineage])

		Doctest:
		>>> test_resovler = TestTaxonNamesResolver(scenario = 1)
		>>> genTaxTree(test_resovler, draw  = False)
		(Tree(weight=1.0, rooted=False), ['claA'])"""
	ranks = resolver.retrieve('classification_path_ranks')
	if by == 'qnames' or by == 'name_string':
		lineages = resolver.retrieve('classification_path')
		if by == 'qnames':
			idents = resolver.retrieve('query_name')
		else:
			idents = resolver.retrieve('name_string')
	elif by == 'taxids':
		idents = resolver.retrieve('taxon_id')
		idents = ["tx{0}".format(e) for e in idents]
		lineages = resolver.retrieve('classification_path_ids')          
	else:
		raise NameError("Unrecognised 'by' argument (either qnames, taxids or name_string)")
	for i, lineage in enumerate(lineages):
		lineage.reverse()
		lineages[i] = lineage
	for i, rank in enumerate(ranks):
		rank.reverse()
		ranks[i] = rank
	# make lineages of same ranks
	all_ranks = [e2 for e1 in ranks for e2 in e1]
	rank_freq =  collections.Counter(all_ranks).items()
	shared_ranks = [e for e, f in rank_freq if f == len(idents)]
	line_bool = [[1 if e2 in shared_ranks else 0 for e2 in e1] for e1 in ranks]
	lineages = [[lineages[i1][i2] for i2, e2 in enumerate(e1) if e2 == 1] for \
					i1, e1 in enumerate(line_bool)]
	all_lines = [e2 for e1 in lineages for e2 in e1]
	line_freq =  collections.Counter(all_lines).items()
	shared_lineage = [e for e, f in line_freq if f == len(idents)]
	# create line_obj, a tuple of ident and lineage
	line_obj = zip(idents, lineages)
	for i in range(len(lineages[0])):
		for uniq in set([each[1][i] for each in line_obj]):
			# find shared taxonomic groups
			new_node = [each[0] for each in line_obj if each[1][i] == uniq]
			if len(new_node) > 1:
				# extract shared lineage
				lineage = [each[1] for each in line_obj if each[0] == new_node[0]]
				# remove shareds from line_obj
				line_obj = [each for each in line_obj if not each[0] in new_node]
				# convert to strings
				new_node = [str(each) for each in new_node]
				new_node = [re.sub("\\s", "_", each) for each in new_node]
				# add new node to line_obj
				new_node = ('(' + ','.join(new_node) + ')' + str(uniq), lineage[0])
				line_obj.append(new_node)
		if len(line_obj) < 1:
			break
	tree = Phylo.read(StringIO(line_obj[0][0] + ';'), "newick")
	if draw:
		Phylo.draw_ascii(tree)
	return tree, shared_lineage

# I'm not sure this works, nor if it's necessary
#def lineageMerge(resolver, min_rank = 'genus'):
#	"""Return TaxonNamesResolver class with duplicated lineages merged
#
#	Arguments:
#	 resovler = TaxonNamesResolver class
#	 min_rank = the minimum rank level below which all resolved names
#	  will be merged (default genus)
#
#	Return:
#	 TaxonNamesResolver class"""
#	def belowMinRank(index):
		# if min_rank is in rank lineage, then name must be below
#		rtest = [e.lower() for e in ranks[i][:-2]] # drop last element
#		return min_rank.lower() in rtest
#	lineage_ids = resolver.retrieve('classification_path_ids')
#	ranks = resolver.retrieve('classification_path_ranks')
#	for i in range(len(lineage_ids)):
#		merge = False
#		references = lineage_ids[:]
#		query = references.pop(i)[-1] # take last element from lineage
#		merge_list = []
#		for j,reference in enumerate(references):
#			if query in reference[:-2]:
#				if not belowMinRank(i): # delete
#					to_delete = resolver.retrieve('query_name')[i]
#					del resolver._store[to_delete]
#					break
#				else: # merge
#					merge = True
#					merge_list.append(j)
#		if merge:
#			namei = resolver.retrieve('query_name')[i]
#			resultsi = resolver._store[namei]
#			for j,to_merge in enumerate(merge_list):
#				resolver._store[to_merge] = resultsi
#	return resolver
			

def extractHighestClade(resolver, by_ids = False):
	"""Return highest unique taxonomic clade for each resolved name.
	
	Arguments:
	 resolver = TaxonNamesResolver class (lineages should be merged)
	 by_ids = return IDs or names (logical, default False)
	 
	Return:
	 (query_name, name or ID)

	Doctest:
	>>> test_resolver = TestTaxonNamesResolver(scenario = 2)
	>>> extractHighestClade(test_resolver, by_ids = True)
	(['term5', 'term4', 'term3', 'term2', 'term1'], [22, 32, 42, 2, 41])"""
	def find(i):
		# for each lineage, step through the lineage
		# keep all reference lineages that share the same lineage
		#  as the query lineage
		# as soon as the query lineage is unique, return its index
		references = lineages[:]
		query = references.pop(i)
		for j,qid in enumerate(query):
			refids = []
			too_small = []
			for k in range(len(references)):
				try:
					refids.append(references[k][j])
				except IndexError: # lineages that are too small must be dropped
					too_small.append(k)
			matches = [e == qid for e in refids]
			references = [e for ei,e in enumerate(references) if ei not in too_small]
			if any(matches) or len(references) == 0:
				references = [references[ei] for ei,e in enumerate(matches) if e]
			else:
				return j
		else:
			return j
	q_names = resolver.retrieve('query_name')
	ranks = resolver.retrieve('classification_path_ranks')
	lineages = resolver.retrieve('classification_path_ids')
	lineages = [[int(e2) for e2 in e1] for e1 in lineages]
	indexes = []
	for i in range(len(lineages)):
		indexes.append(find(i))
	if not by_ids:
		lineages = resolver.retrieve('classification_path')
	res = [lineages[ei][e] for ei,e in enumerate(indexes)]
	return q_names, res

if __name__ == '__main__':
	import doctest
	doctest.testmod()
