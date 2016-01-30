#! /usr/bin/env python
# D.J. Bennett
# 24/03/2014
"""
misc tools
"""
from __future__ import absolute_import

# PACKAGES
import re
from six.moves import range

# GLOBALS
# nodes on a taxonomic tree
default_taxonomy = ['subspecies', 'species', 'subgenus', 'genus', 'tribe',
                    'subfamily', 'family', 'superfamily', 'parvorder',
                    'infraorder', 'suborder', 'order', 'superorder',
                    'parvclass', 'infraclass', 'subclass', 'class',
                    'superclass', 'subphylum', 'phylum', 'kingdom',
                    'superkingdom']


# CLASSES
class TaxRef(object):
    '''Reference for taxonimic identities'''
    def __init__(self, ident, rank, taxonomy=default_taxonomy):
        super(TaxRef, self).__setattr__('taxonomy', taxonomy)
        super(TaxRef, self).__setattr__('ident', ident)
        super(TaxRef, self).__setattr__('rank', rank)
        # level is a numerical value for rank in taxonomy
        super(TaxRef, self).__setattr__('level',
                                        self._getLevel(rank, taxonomy))
        super(TaxRef, self).__setattr__('counter', 0)  # count ident changes

    def change(self, ident, rank=None):
        '''Change ident'''
        self.ident = ident
        if rank:
            self.rank = rank
            self.level = self._getLevel(rank, self.taxonomy)
        # count changes made to instance
        self.counter += 1

    def _getLevel(self, rank, taxonomy):
        if rank in taxonomy:
            return taxonomy.index(rank)
        # else find its closest by using the default taxonomy
        dlevel = default_taxonomy.index(rank)
        i = 1
        d = dlevel + i
        up = True
        while i <= len(default_taxonomy):
            if d > 0:
                try:
                    drank = default_taxonomy[d]
                except IndexError:
                    pass
                if drank in taxonomy:
                    return taxonomy.index(drank)
            if up:
                d = dlevel - i
                up = False
            else:
                i += 1
                d = dlevel + i
                up = True

    def __repr__(self):
        return self.ident

    def __str__(self):
        return '{0} -- {1}@{2}'.format(self.ident, self.rank, self.level)

    def __setattr__(self, name, value):
        if name in ['ident', 'rank']:
            if not isinstance(value, str):
                raise ValueError('[{0}] must be a string'.format(name))
            super(TaxRef, self).__setattr__(name, value)
        elif name in ['level', 'counter']:
            if not isinstance(value, int):
                raise ValueError('[{0}] must be an integer'.format(name))
            super(TaxRef, self).__setattr__(name, value)
        else:
            raise AttributeError(name)


class TaxDict(dict):
    '''Taxonomic Dictionary : hold and return taxonomic information'''

    def __init__(self, idents, ranks, lineages, taxonomy=default_taxonomy,
                 **kwargs):
        # add entry for each ident of lineages ordered by taxonomy
        # ranks without corresponding lineage are given ''
        # 'ident' is the unique name for a taxonomic entity (e.g. query name)
        # 'ranks' must be the names of the corresponding ranks in lineages
        #  (e.g. classification_path_ranks)
        # 'lineages' is the names for each of the ranks (e.g.
        #   classification_path or classification_path_ids)
        if taxonomy:
            self.taxonomy = taxonomy
        else:
            self.taxonomy = default_taxonomy
        for i in range(len(idents)):
            # extract lineage according to given taxonomy
            lineage = [lineages[i][ranks[i].index(e)] if e in ranks[i] else ''
                       for e in self.taxonomy]
            # create taxref
            taxref = TaxRef(ident=idents[i], rank=ranks[i][-1],
                            taxonomy=self.taxonomy)
            # create key for ident and insert a dictionary of:
            #  lineage, taxref, cident, ident and rank
            self[idents[i]] = {'lineage': lineage, 'taxref': taxref,
                               'cident': None, 'rank': ranks[i][-1],
                               'ident': lineage[taxref.level]}
        # add addtional optional slots from **kwargs
        self._additional(idents, kwargs)
        # gen hierarchy
        self._hierarchy()
        # contexualise
        self._contextualise()

    def _additional(self, idents, kwargs):
        '''Add additional data slots from **kwargs'''
        if kwargs:
            for name, value in list(kwargs.items()):
                if not isinstance(value, list):
                    raise ValueError('Additional arguments must be lists of \
same length as idents')
                for i in range(len(value)):
                    self[idents[i]][name] = value[i]

    def _slice(self, level):
        '''Return list of tuples of ident and lineage ident for given level
(numbered rank)'''
        if level >= len(self.taxonomy):
            raise IndexError('Level greater than size of taxonomy')
        res = []
        for ident in sorted(list(self.keys())):
            res.append((self[ident]['taxref'], self[ident]['lineage'][level]))
        return res

    def _group(self, taxslice):
        '''Return list of lists of idents grouped by shared rank'''
        res = []
        while taxslice:
            taxref, lident = taxslice.pop()
            if lident == '':
                res.append(([taxref], lident))
            else:
                # identify idents in the same group and pop from taxslice
                i = 0
                group = []
                while i < len(taxslice):
                    if taxslice[i][1] == lident:
                        group.append(taxslice.pop(i)[0])
                    else:
                        i += 1
                group.append(taxref)
                res.append((group, lident))
        return res

    def _hierarchy(self):
        '''Generate dictionary of referenced idents grouped by shared rank'''
        self.hierarchy = {}
        for rank in self.taxonomy:
            # extract lineage idents for this rank
            taxslice = self._slice(level=self.taxonomy.index(rank))
            # group idents by shared group at this rank
            self.hierarchy[rank] = self._group(taxslice)

    def _contextualise(self):
        '''Determine contextual idents (cidents)'''
        # loop through hierarchy identifying unique lineages
        # TODO: gain other contextual information, not just ident
        deja_vues = []
        for rank in reversed(self.taxonomy):
            # return named clades -- '' are ignored
            clades = [e for e in self.hierarchy[rank] if e[1]]
            # print 'Rank: {0} - {1}'.format(rank, len(clades))
            # get unique lineages at this level
            uniques = [e for e in clades if len(e[0]) == 1]
            # removed those already seen
            uniques = [e for e in uniques if e[0][0].ident not in deja_vues]
            # add each to self[ident]['cident']
            for e in uniques:
                ident = e[0][0].ident
                self[ident]['cident'] = e[1]
                deja_vues.append(ident)


# FUNCTIONS
def stringClade(taxrefs, name, at):
    '''Return a Newick string from a list of TaxRefs'''
    string = []
    for ref in taxrefs:
        # distance is the difference between the taxonomic level of the ref
        #  and the current level of the tree growth
        d = float(at-ref.level)
        # ensure no spaces in ident, Newick tree cannot have spaces
        ident = re.sub("\s", "_", ref.ident)
        string.append('{0}:{1}'.format(ident, d))
    # join into single string with a name for the clade
    string = ','.join(string)
    string = '({0}){1}'.format(string, name)
    return string


def taxTree(taxdict):
    """Return taxonomic Newick tree"""
    # the taxonomic dictionary holds the lineage of each ident in
    #  the same order as the taxonomy
    # use hierarchy to construct a taxonomic tree
    for rank in taxdict.taxonomy:
        current_level = float(taxdict.taxonomy.index(rank))
        # get clades at this rank in hierarchy
        clades = taxdict.hierarchy[rank]
        # merge those that are in the same clade into a cladestring
        for clade in clades:
            # unpack the identities in this clade and its clade name
            cladeidents, cladename = clade
            # Remove '' TaxRefs -- in cladestring already
            cladeidents = [e for e in cladeidents if e.ident]
            # only create cladestring if more than one ident in clade
            if len(cladeidents) < 2:
                continue
            # label node by 'clade'_'rank'
            cladename = '{0}_{1}'.format(cladename, rank)
            cladestring = stringClade(cladeidents, cladename, current_level)
            # replace first TaxRef in cladeidents with cladestring
            cladeidents[0].change(ident=cladestring, rank=rank)
            # replace all other TaxRefs with ''
            for e in cladeidents[1:]:
                e.change(ident='', rank=rank)
    # join any remaining strands into tree
    if len(taxdict.hierarchy[taxdict.taxonomy[-1]]) > 1:
        # unlist first
        clade = [e[0] for e in taxdict.hierarchy[taxdict.taxonomy[-1]]]
        cladeidents = sum(clade, [])
        cladeidents = [e for e in cladeidents if e.ident]
        cladestring = stringClade(cladeidents, 'life', current_level+1)
    return cladestring + ';'
