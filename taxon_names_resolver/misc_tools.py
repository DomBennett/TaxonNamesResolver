#! /usr/bin/env python
# D.J. Bennett
# 24/03/2014
"""
pglt names tools
"""

# PACKAGES
import re

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

    def __init__(self, idents, ranks, lineages, taxonomy=default_taxonomy):
        # add entry for each ident of lineages ordered by taxonomy
        # ranks without corresponding lineage are given ''
        self.taxonomy = taxonomy
        for i in range(len(idents)):
            # extract lineage according to given taxonomy
            lineage = [lineages[i][ranks[i].index(e)] if e in ranks[i] else ''
                       for e in taxonomy]
            # create taxref
            taxref = TaxRef(ident=idents[i], rank=ranks[i][-1],
                            taxonomy=taxonomy)
            # create key for ident and insert dictionary of lineage and taxref
            self[idents[i]] = {'lineage': lineage, 'taxref': taxref}

    def _slice(self, level):
        '''Return list of tuples of ident and lineage ident for given level
(numbered rank)'''
        if level >= len(self.taxonomy):
            raise IndexError('Level greater than size of taxonomy')
        res = []
        for ident in self.keys():
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

    def hierarchy(self):
        '''Return dictionary of referenced idents grouped by shared rank'''
        res = {}
        for rank in self.taxonomy:
            # extract lineage idents for this rank
            taxslice = self._slice(level=self.taxonomy.index(rank))
            # group idents by shared group at this rank
            res[rank] = self._group(taxslice)
        return res


# FUNCTIONS
def stringClade(taxrefs, name, at):
    '''Return a Newick string from a list of TaxRefs'''
    string = []
    for ref in taxrefs:
        # distance is the difference between the taxonomic level of the ref
        #  and the current level of the tree growth
        d = float(at-ref.level)
        string.append('{0}:{1}'.format(ref.ident, d))
    # join into single string with a name for the clade
    string = ','.join(string)
    string = '({0}){1}'.format(string, name)
    return string


def taxTree(idents, ranks, lineages, taxonomy=None):
    """Generate Taxonomic Newick tree"""
    if not taxonomy:
        taxonomy = default_taxonomy
    # replace any ' ' with '_' for taxon tree
    idents = [re.sub("\s", "_", e) for e in idents]
    # create a taxonomic dictionary which holds the lineage of each ident in
    #  the same order as the given taxonomy
    taxdict = TaxDict(idents=idents, ranks=ranks, lineages=lineages,
                      taxonomy=taxonomy)
    # generate a hierarchy from the taxdict
    hierarchy = taxdict.hierarchy()
    # use hierarchy to construct a taxonomic tree
    for rank in taxonomy[1:]:
        current_level = float(taxonomy.index(rank))
        # get clades at this rank in hierarchy
        clades = hierarchy[rank]
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
    if len(hierarchy[taxonomy[-1]]) > 1:
        # unlist first
        clade = [e[0] for e in hierarchy[taxonomy[-1]]]
        cladeidents = sum(clade, [])
        cladeidents = [e for e in cladeidents if e.ident]
        cladestring = stringClade(cladeidents, 'life', current_level+1)
    return cladestring + ';'
