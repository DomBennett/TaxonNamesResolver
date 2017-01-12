#! /usr/bin/env python
# D. J. Bennett
# 16/05/2014
"""
Tools for interacting with the GNR.
"""
from __future__ import absolute_import

import time
import contextlib
import json
import os
import six
from six.moves import urllib


# FUNCTIONS
def safeReadJSON(url, logger, max_check=6, waittime=30):
    '''Return JSON object from URL'''
    counter = 0
    # try, try and try again ....
    while counter < max_check:
        try:
            with contextlib.closing(urllib.request.urlopen(url)) as f:
                res = json.loads(f.read().decode('utf8'))
            return res
        except Exception as errmsg:
            logger.info('----- GNR error [{0}] : retrying ----'.format(errmsg))
        counter += 1
        time.sleep(waittime)
    logger.error('----- Returning nothing : GNR server may be down -----')
    return None


# CLASSES
class GnrDataSources(object):
    """GNR data sources class: extract IDs for specified data sources."""

    def __init__(self, logger):
        url = 'http://resolver.globalnames.org/data_sources.json'
        self.available = safeReadJSON(url, logger)

    def summary(self):
        # see what sources are available
        return [dict(id=ds['id'], title=ds['title']) for ds in self.available]

    def byName(self, names, invert=False):
        if invert:
            return [ds['id'] for ds in self.available if not ds['title'] in
                    names]
        else:
            return [ds['id'] for ds in self.available if ds['title'] in names]


class GnrResolver(object):
    """GNR resolver class: search the GNR"""

    def __init__(self, logger, datasource='NCBI'):
        self.logger = logger
        ds = GnrDataSources(logger)
        self.write_counter = 1
        self.Id = ds.byName(datasource)
        self.otherIds = ds.byName(datasource, invert=True)
        self.waittime = 600  # wait ten minutes if server fail
        self.max_check = 6  # search for up to an hour

    def search(self, terms, prelim=True):
        """Search terms against GNR. If prelim = False, search other datasources \
for alternative names (i.e. synonyms) with which to search main datasource.\
Return JSON object."""
        # TODO: There are now lots of additional data sources, make additional
        # searching optional (11/01/2017)
        if prelim:  # preliminary search
            res = self._resolve(terms, self.Id)
            self._write(res)
            return res
        else:  # search other DSs for alt names, search DS with these
            # quick fix: https://github.com/DomBennett/TaxonNamesResolver/issues/5
            # seems to be due to limit on number of ids in single request
            # switiching to a for loop for each data source
            # appending all results into single res
            res = []
            for ds_id in self.otherIds:
                tmp = self._resolve(terms, [ds_id])
                res.append(tmp[0])
            self._write(res)
            alt_terms = self._parseNames(res)
            if len(alt_terms) == 0:
                return False
            else:
                # search the main source again with alt_terms
                # replace names in json
                terms = [each[1] for each in alt_terms]  # unzip
                res = self._resolve(terms, self.Id)
                self._write(res)
                alt_res = self._replaceSupStrNames(res, alt_terms)
                return alt_res

    def _parseNames(self, jobj):
        # return a list of tuples (term, name) from second search
        # TODO(07/06/2013): record DSs used
        alt_terms = []
        for record in jobj:
            if 'results' not in list(record.keys()):
                pass
            else:
                term = record['supplied_name_string']
                results = record['results']
                for result in results:
                    r_name = result['canonical_form']
                    if r_name == term or r_name is None:
                        continue
                    alt_terms.append((term, r_name))
        alt_terms = list(set(alt_terms))
        return alt_terms

    def _replaceSupStrNames(self, jobj, alt_terms):
        # replace sup name in jobj with original terms
        for record in jobj:
            sup_name = record['supplied_name_string']
            # find original name in alt_terms
            term = [i for i, each in enumerate(alt_terms) if each[1] ==
                    sup_name]
            # pop from alt_terms and rename json use 0 to
            #  avoid the possibility of having the same term with >1 r_names
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
            self.logger.info('Querying [{0}] to [{1}] of [{2}]'.
                             format(lower, upper, len(terms)))
            query = self._query(terms[lower:upper], ds_id)
            res.append(query)
            lower = upper
        res = [record for search in res for record in search['data']]
        return(res)

    def _query(self, terms, data_source_ids):
        ds_ids = [str(id) for id in data_source_ids]
        terms = [urllib.parse.quote(six.text_type(t).encode('utf8')) for t in terms]
        url = ('http://resolver.globalnames.org/name_resolvers.json?' +
               'data_source_ids=' + '|'.join(ds_ids) + '&' +
               'resolve_once=false&' + 'names=' + '|'.join(terms))
        return safeReadJSON(url, self.logger)

    def _write(self, jobj):
        directory = os.path.join(os.getcwd(), 'resolved_names')
        filename = "{0}_raw_results.json".format(self.write_counter)
        jobj_file = os.path.join(directory, filename)
        with open(jobj_file, 'w') as outfile:
            json.dump(jobj, outfile)
        self.write_counter += 1


class GnrStore(dict):
    """GNR store class: acts like a dictionary for GNR JSON format"""

    def __init__(self, terms, logger, tax_group=None):
        self.logger = logger
        # Issue 6: suggest multiple tax_groups, not just one
        if not tax_group:
            self.tax_group = tax_group
        else:
            if not isinstance(tax_group, list):
                tax_group = [tax_group]
            # ensure strings
            self.tax_group = [str(e) for e in tax_group]
        for term in terms:
            self[term] = []

    def _filter(self, results):
        # filter out all results that are not in tax_group
        if not self.tax_group:
            return results
        filtered = []
        for result in results:
            classids = result['classification_path_ids'].split('|')
            if any([True if e in classids else False for e in self.tax_group]):
                filtered.append(result)
        return filtered

    def add(self, jobj):
        if not isinstance(jobj, bool):
            for record in jobj:
                term = record['supplied_name_string']
                try:
                    if 'results' in list(record.keys()):
                        results = self._filter(record['results'])
                        self[term].extend(results)
                except KeyError:
                    self.logger.debug('JSON object contains terms not in self.logger')

    def replace(self, jobj):
        for record in jobj:
            term = record['supplied_name_string']
            try:
                if 'results' in list(record.keys()):
                    results = self._filter(record['results'])
                    self[term] = results
                else:
                    self[term] = []
            except KeyError:
                self.logger.debug('JSON object contains terms not in GnrStore')
