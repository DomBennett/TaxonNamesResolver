#!/usr/bin/python
## No warranty, no copyright
## Dominic John Bennett
## 03/12/2013
## TODO: adapt genTaxTree to use taxonomic information outside of named ranks
## TODO: write doctests
## TODO: test mergeLineage and extractHighestClade with duplicated test data

import collections, re
from Bio import Phylo
from cStringIO import StringIO

def genTaxTree(resolver, by = 'qnames', draw = False):
	"""Generate Newick tree from TaxonNamesResolver class.
        
        Arguments:
         resolver = TaxonNamesResolver class
         by = Tip labels, either 'qnames', 'taxids' or 'name_string'
         draw = Draw ascii tree (logical)

         Return:
          (Newick Tree Object, [shared lineage])"""
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
		raise Error("Unrecognised 'by' argument (either qnames, taxids or name_string)")
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

def lineageMerge(resolver, min_rank = 'genus'):
    """Return TaxonNamesResolver class with duplicated lineages merged

    Arguments:
     resovler = TaxonNamesResolver class
     min_rank = the minimum rank level below which all resolved names
      will be merged (default genus)

    Return:
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
            

def extractHighestClade(resolver, by_ids = False):
    """Return highest unique taxonomic clade for each resolved name.
    
    Arguments:
     resolver = TaxonNamesResolver class (lineages should be merged)
     by_ids = return IDs or names (logical, default False)
     
    Return:
     (query_name, name or ID)"""
    def find(i):
        # for each lineage, step through the lineage
        # keep all reference lineages that share the same lineage
        #  as the query lineage
        # as soon as the query lineage is unique, return its index
        references = lineages[:]
        query = references.pop(i)
        for j,qid in enumerate(query):
            refids = [e[j] for e in references]
            matches = [e == qid for e in refids]
            if any(matches):
                references = [references[ei] for ei,e in enumerate(matches) if e]
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
	pass # doctests to come
