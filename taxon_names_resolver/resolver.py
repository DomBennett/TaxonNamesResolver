#! /usr/bin/env python
# D.J. Bennett
# 16/05/2014
# TODO: use a minimal score (avoid calling hybrids 'species')
# TODO: create a stats-out function
"""
Resolver class for parsing GNR records.
"""
from __future__ import absolute_import

import json
import os
import csv
import re
import copy
import logging
from .gnr_tools import GnrStore
from .gnr_tools import GnrResolver
import six
from six.moves import zip, urllib


# CLASSES
class EncodingError(Exception):
    pass


class Resolver(object):
    """Taxon Names Resovler class : Automatically resolves taxon names \
through GNR. All output written in 'resolved_names' folder.
See https://github.com/DomBennett/TaxonNamesResolver for details."""

    def __init__(self, input_file=None, datasource='NCBI', taxon_id=None,
                 terms=None, lowrank=False, logger=logging.getLogger('')):
        # add logger
        self.logger = logger
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
                self.logger.info("No terms provided")
        terms = list(set(terms))
        self.logger.info('Found [{0}] taxon names to search in input file... '.
                         format(len(terms)))
        self.logger.info('... of which [{0}] are unique.'.format(len(terms)))
        # init dep classes
        self._check(terms)
        self.terms = terms
        self._res = GnrResolver(logger=self.logger, datasource=datasource)
        self.primary_datasource = datasource
        self._store = GnrStore(terms, tax_group=taxon_id, logger=self.logger)
        # http://resolver.globalnames.org/api
        self.key_terms = ['query_name', 'classification_path',
                          'data_source_title', 'match_type', 'score',
                          'classification_path_ranks', 'name_string',
                          'canonical_form', 'classification_path_ids',
                          'prescore', 'data_source_id', 'taxon_id', 'gni_uuid']
        self.lowrank = lowrank  # return lowest ranked match
        # self.tnr_obj = [] # this will hold all output

    def _check(self, terms):
        """Check terms do not contain unknown characters"""
        for t in terms:
            try:
                _ = urllib.parse.quote(six.text_type(t).encode('utf8'))
            except:
                self.logger.error('Unknown character in [{0}]!'.format(t))
                self.logger.error('.... remove character and try again.')
                raise EncodingError

    def main(self):
        """Search and sieve query names."""
        # TODO: Break up, too complex
        primary_bool = True
        no_records = True
        nsearch = 1
        search_terms = self.terms
        original_names = []
        while True:
            if primary_bool:
                self.logger.info('Searching [{0}] ...'.format(
                    self.primary_datasource))
            else:
                self.logger.info('Searching other datasources ...')
            res = self._res.search(search_terms, prelim=primary_bool)
            if nsearch > 2 and res:
                # if second search failed, look up alternative names
                for each_res, original_name in zip(res, original_names):
                    each_res['supplied_name_string'] = original_name
            self._store.add(res)
            # Check for returns without records
            no_records = self._count(nrecords=1)
            if no_records:
                if nsearch == 1:
                    primary_bool = False
                elif nsearch == 2:
                    original_names = no_records
                    # genus names
                    no_records = [e.split()[0] for e in no_records]
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
        multi_records = self._count(greater=True, nrecords=1)
        if multi_records:
            self.logger.info('Choosing best records to return ...')
            res = self._sieve(multi_records)
            self._store.replace(res)

    # def extract(self, what): # depends on tnr
    #    lkeys = ['qnames', 'rnames', 'taxonids', 'ranks']
    #    i = [i for i, each in enumerate(lkeys) if what is each][0]
    #    res = [each[i] for each in self.tnr_obj]
    #    return res

    def _readInJson(self):
        jobj_file = os.path.join(self.outdir, 'prelim_search_results.json')
        with open(jobj_file, 'r') as infile:
            jobj = json.load(infile)
        return jobj

    def _count(self, greater=False, nrecords=0):
        # return a list of all keys that have records
        #  greater or less than nrecords
        GnrStore = self._store
        assessed = []
        lens = [len(GnrStore[key]) for key in list(GnrStore.keys())]
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
interest, have the highest score, have the lowest taxonomic rank (if lowrank is
true) and/or are the first item returned."""
        # TODO: Break up, too complex
        GnrStore = self._store

        def writeAsJson(term, results):
            record = {'supplied_name_string': term}
            if len(results) > 0:
                record['results'] = results
            return record

        def boolResults(results, bool_li, rand=False):
            if rand:
                # choose first record (most likely best?)
                results = [results[0]]
            elif sum(bool_li) == 1:
                results = [results[bool_li.index(1)]]
            elif sum(bool_li) == 0:
                # return 'no_record'
                return []
            else:
                results = [result for i, result in enumerate(results) if
                           bool_li[i]]
            return results

        sieved = []
        ranks = ['species', 'genus', 'family', 'order', 'superorder', 'class',
                 'superclass', 'subphylum', 'phylum', 'kingdom',
                 'superkingdom']
        for term in multiple_records:
            results = GnrStore[term]
            while len(results) > 1:
                # choose result with best score
                scores = [result['score'] for result in results]
                bool_score = [1 if score == max(scores) else 0 for score in
                              scores]
                results = boolResults(results, bool_score)
                # choose result resolved to lowest taxonomic rank
                if self.lowrank:
                    res_ranks = [result['classification_path_ranks'].
                                 split('|') for result in results]
                    # calculate 'rank scores' for named and un-named ranks
                    nmd_rnks = []
                    unnmd_rnks = []
                    for rs in res_ranks:
                        nmd_rnks.append(min([j for j,e in enumerate(ranks) if
                                             e in rs]))
                        unnmd_rnk = [j for j,e in enumerate(rs) if
                                     e == ranks[nmd_rnks[-1]]][0]
                        unnmd_rnk -= len(rs)
                        unnmd_rnks.append(unnmd_rnk)
                    # calculate bool
                    unnmd_rnks = [e if nmd_rnks[j] == min(nmd_rnks) else 0 for
                                  j,e in enumerate(unnmd_rnks)]
                    bool_rank = [1 if e == min(unnmd_rnks) else 0 for e in
                                 unnmd_rnks]
                    results = boolResults(results, bool_rank)
                results = boolResults(results, [], rand=True)
            record = writeAsJson(term, results)
            sieved.append(record)
        return sieved

    def write(self):
        """Write csv file of resolved names and txt file of unresolved names.
        """
        csv_file = os.path.join(self.outdir, 'search_results.csv')
        txt_file = os.path.join(self.outdir, 'unresolved.txt')
        headers = self.key_terms
        unresolved = []
        with open(csv_file, 'w') as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            for key in list(self._store.keys()):
                results = self._store[key]
                if len(results) == 0:
                    unresolved.append(key)
                else:
                    row = [key]
                    for key_term in headers[1:]:
                        element = results[0][key_term]
                        # GNR returns UTF-8, csv requires ascii
                        #
                        # *** Note ***
                        # According to all docs for csv versions >= 2.6, csv
                        # can handle either UTF-8 or ascii, just not Unicode.
                        # In py3, the following two lines result in csv printing
                        # the element with a bitstring. If GNR is actually
                        # returning UTF-8, it seems easiest to just drop these

                        # if 'encode' in dir(element):
                        #     element = element.encode('ascii')
                        row.append(element)
                    writer.writerow(row)
        if len(unresolved) > 0:
            with open(txt_file, 'w') as file:
                for name in unresolved:
                    file.write("{0}\n".format(name))

    def retrieve(self, key_term):
        """Return data for key term specified for each resolved name as a list.
Possible terms (02/12/2013): 'query_name', 'classification_path',
'data_source_title', 'match_type', 'score', 'classification_path_ranks',
'name_string', 'canonical_form',\
'classification_path_ids', 'prescore', 'data_source_id', 'taxon_id',
'gni_uuid'"""
        if key_term not in self.key_terms:
            raise IndexError('Term given is invalid! Check doc string for \
valid terms.')
        store = self._store
        retrieved = []
        for key in list(store.keys()):
            # take copy, so changes made to the returned list do not affect
            #  store
            record = copy.deepcopy(store[key])
            if len(record) > 0:
                if key_term == 'query_name':
                    retrieved.append(key)
                else:
                    retrieved.append(record[0][key_term])
        if re.search('path', key_term):
            retrieved = [[r2 for r2 in r1.split('|')[1:]] for r1 in retrieved]
        return retrieved
