# TaxonNamesResolver
Automatically search taxon names against the [Global Names Resolver (GNR)](resolver.globalnames.org) to:
* Find the most likely matching name for a given taxonomic name
* Generate taxonomic lineages
* Generate taxonomic IDs

Return raw JSON search files, csv file of resolved names and txt file of unresolved names.

## How does it work?
Global Names Resolver is a powerful API that uses specialised fuzzy-matching algorithms to search taxonomic datasources across the web. TaxonNamesResovler queries these datasources through the programming language python.

TaxonNamesResolver first searches all the names (in chunks of 100) against the main datasource (search 1). Names that fail to be resolved are then searched against other datasources (search 2) to find synonyms. If any synonyms are returned, these are searched against the main datasource (search 3). Names that remain unresolved, are reduced to their genus name, and these are again searched as above (searches 4 to 6).

For returned names with multiple records, the most likely match is found by testing if the name is in the correct clade (as specified by the user), if the name has the highest GNR score and/or the name at the lowest taxonomic level.

The resolved names are then written to a csv file.

## Basic usage
Requires [python 2.7.3](https://wiki.python.org/moin/BeginnersGuide/Download) to be installed and the installation of the following python packages: contextlib, json, urllib and urllib2
* Download taxon_names_resolver.py
* Place in folder with list of taxon names to be searched
* Open terminal/command prompt and type: 'python taxon_names_resolver.py'
 * Provide file name for list of taxon names
 * Provide the name of the datasource you wish to search against (NCBI by default)
 * Provide the taxon ID of the parent clade of all names given (see Taxon ID below)
 * Ensure details are correct, then click enter
* Results are generated in the 'resolved_names' folder

## Within python
The TaxonNamesResolver class can be imported into a python session and run, so:
```{python}
resolver = TaxonNamesResolver(input_file, datasource, taxon_id)
resolver.main() # to run the search
resolver.write() # to output the csv file
```
Different elements of the returned JSON file can be extracted using `.retrieve(key_term)`. This function will return a list for each resolved name based on the key_term given. The key_terms are:
* query_name -- name used to search datasources
* classification_path -- taxonomic lineage, according to main datasource
* data_source_title -- name of datasource
* match_type -- see http://resolver.globalnames.org/api
* score -- see http://resolver.globalnames.org/api
* classification_path_ranks -- the rank name for the lineage, no name specified by ''
* name_string -- the returned name in the datasorce
* canonical_form -- the name's canonical form (i.e. Genus species)
* classification_path_ids -- the taxonomic IDs for the lineage (datasource dependent)
* prescore -- see http://resolver.globalnames.org/api
* data_source_id -- GNR ID for datasource
* taxon_id -- taxonomic ID for the returned name (datasource dependent)
* gni_uuid -- Global Names Unique ID

## Taxon ID
Taxonomy is full of synonyms; to avoid returning the wrong name it is best to specifiy a parent clade ID. For example, if you know all the taxon names in your list are mammals, then any names matched that are not a mammal must be incorrect and TaxonNamesResolver will remove them. When searching NCBI as the main datasource, use [NCBI taxonomy](http://www.ncbi.nlm.nih.gov/taxonomy) to search for and find the ID of your parent clade.

## JSON files
The GNR API returns [JSON](http://en.wikipedia.org/wiki/JSON) files after searching across datasources. To keep the program transparent, every search carried out be TaxonNamesResolver is saved in the 'resolved_names' folder. These files can be viewed in a web browser with the appropriate add-on.

## Future
* Unit tests
* TaxonNamesResolver tools:
 * Generate taxonomic trees
 * Return different names based on given options

## Thanks
* Thanks to Lawrence Hudson for coding inspiration
* Thanks to the GNR API and the [Global Names](http://www.globalnames.org/) family for creating this wonderful resource

## Author
Dominic John Bennett (dominic.john.bennett@gmail.com)