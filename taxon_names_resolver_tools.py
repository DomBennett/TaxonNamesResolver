#!/usr/bin/python
## No warranty, no copyright
## Dominic John Bennett
## 03/12/2013
## TODO: adapt genTaxTree to use taxonomic information outside of named ranks
## TODO: write doctests
## TODO: test mergeLineage and extractHighestClade with duplicated test data

import collections
from Bio import Phylo
from cStringIO import StringIO

def genTaxTree(resolver, by_ids = False, draw = False):
	"""Generate Newick tree from TaxonNamesResolver class.
        
        Arguments:
         resolver = TaxonNamesResolver class
         by_ids = Use taxon IDs instead of names (logical)
         draw = Draw ascii tree (logical)

         Returns:
          (Newick Tree Object, [shared lineage])"""
        if by_ids:
            idents = resolver.retrieve('taxon_id')
            lineages = resolver.retrieve('classification_path_ids')
            ranks = resolver.retrieve('classification_path_ranks')
        else:
            idents = resolver.retrieve('name_string')
            lineages = resolver.retrieve('classification_path')
            ranks = resolver.retrieve('classification_path_ranks')
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
	lineages = [[lineages[i1][i2] for i2, e2 in enumerate(e1) if e2 == 1] for i1, e1 in enumerate(line_bool)]
	all_lines = [e2 for e1 in lineages for e2 in e1]
	line_freq =  collections.Counter(all_lines).items()
	shared_lineage = [e for e, f in line_freq if f == len(idents)]
	# TODO: if shared lineage is empty... drop radically different taxa
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
				# add new node to line_obj
				new_node = ('(' + ','.join(new_node) + ')', lineage[0])
				line_obj.append(new_node)
		if len(line_obj) < 1:
			break
        tree = Phylo.read(StringIO(line_obj[0][0] + ';'), "newick")
        if draw:
            Phylo.draw_ascii(tree)
	return (tree, shared_lineage)

def lineageMerge(resovler, min_rank = 'genus'):
    """Return TaxonNamesResolver class with duplicated lineages merged

    Arguments:
     resovler = TaxonNamesResolver class
     min_rank = the minimum rank level below which all resolved names
      will be merged (default genus)

    Returns:
     TaxonNamesResolver class"""
    def belowMinRank(index):
        # if min_rank is in rank lineage, then name must be below
        rtest = [e.lower() for e in ranks[i][:-2]] # drop last element
        return min_rank.lower() in rtest
    lineage_ids = resolver.retrieve('classification_path_ids')
    ranks = resolver.retrieve('classification_path_ranks')
    for i in range(len(lineage_ids)):
        merge = False
        references = lineage_ids[:]
        query = references.pop(i)[-1] # take last element from lineage
        merge_list = []
        for j,reference in enumerate(references):
            if query in reference[:-2]:
                if not belowMinRank(i): # delete
                    to_delete = resolver.retrieve('query_name')[i]
                    del resolver._store[to_delete]
                    break
                else: # merge
                    merge = True
                    merge_list.append(j)
        if merge:
            namei = resolver.retrieve('query_name')[i]
            resultsi = resolver._store[namei]
            for j,to_merge in enumerate(merge_list):
                resolver._store[to_merge] = resultsi
    return resolver
            

def extractHighestGroup(resolver):
    """Return highest unique taxonomic group for each resolved name."""
    def find(i):
        temp = lineages[:]
        lineage = temp.pop(i)
        for j, lid in enumerate(lineage):
            tempslice = [e[j] for e in temp]
            matches = [e == lid for e in tempslice]
            if any(matches):
                temp = [temp[ei] for ei,e in enumerate(matches) if e]
            else:
                return j,lid
            q_names = self.retrieve('query_name')
            lineages_id = self.retrieve('classification_path_ids')
            lineages_id = [int(e) for e in lineages_id]
            lineages = self.retrieve('classification_path')
            ranks = self.retrieve('classification_path_ranks')
            ## TODO: remove duplicated resovled taxa
            for i in len(lineages):
                j,lid = find(i)
            self.lineage_ids = lineage_ids
            self.huids = huids
            self.ranks = ranks

if __name__ == '__main__':    
    # doctests to come
    merged_resolver = lineageMerge(resolver)
