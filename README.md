# TaxonNamesResolver
Automatically search taxon names against the [Global Names Resolver (GNR)](resolver.globalnames.org) to find the most likely matching name for a given taxonomic name, generate taxonomic lineages and taxonomic IDs.

Returns raw JSON search files, .csv files of resolved names and .txt files of unresolved names.

## Installation
Requires [python 2.7](https://wiki.python.org/moin/BeginnersGuide/Download) to be installed and requires the installation of the python package `setuptools`.
* Click 'Download ZIP' to download the program.
* Open a terminal/command prompt, and navigate to the downloaded folder.
* Type `python setup.py install` (for Unix systems you may need to type `sudo python setup.py install`).
* This will install taxon names resolver as a package in your python library, it will also add the python script `TaxonNamesResolver.py` to your system path.

## How does it work?
Global Names Resolver is a powerful API that uses specialised fuzzy-matching algorithms to search taxonomic datasources across the web. TaxonNamesResovler queries these datasources through the programming language python.

TaxonNamesResolver first searches all the names (in chunks of 100) against the main datasource (search 1). Names that fail to be resolved are then searched against other datasources (search 2) to find synonyms. If any synonyms are returned, these are searched against the main datasource (search 3). Names that remain unresolved, are reduced to their genus name, and these are again searched as above (searches 4 to 6).

For returned names with multiple records, the most likely match is found by testing if the name is in the correct clade (as specified by the user), if the name has the highest GNR score and/or the name at the lowest taxonomic level.

The resolved names are then written to a .csv file.

### Taxon ID
Taxonomy is full of synonyms; to avoid returning the wrong name it is best to specify a parent clade ID. For example, if you know all the taxon names in your list are mammals, then any names matched that are not a mammal must be incorrect and TaxonNamesResolver will remove them. When searching NCBI as the main datasource, use [NCBI taxonomy](http://www.ncbi.nlm.nih.gov/taxonomy) to search for and find the ID of your parent clade. (For mammals it would be 40674.)

### JSON files
The GNR API returns [JSON](http://en.wikipedia.org/wiki/JSON) files after searching across datasources. To keep the program transparent, every search carried out be TaxonNamesResolver is saved in the 'resolved_names' folder. These files can be viewed in a web browser with the appropriate add-on.

## Quick guide
### Command-line usage
Open a terminal/command prompt and type: `TaxonNamesResolver.py`. The program will then ask for: a filepath to a list of names, the datasource you wish to search against and a taxonomic ID of a parent clade for names given (optional). It will then search online and write results to 'resolved_names' folder in the working directory.

### Within python
#### Running
The `taxon_names_resolver` is a python package and can be imported into a session and run, so:
```{python}
from taxon_names_resolver import Resolver
resolver = Resolver(input_file, datasource, taxon_id)
resolver.main()  # to run the search
resolver.write()  # to output the csv file
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

#### Generating a taxonomic tree
`taxTree` creates a Newick string from lists of names, lineages and ranks. It can be used in conjunction with `Resolver`, like:

```{python}
# IMPORT
from taxon_names_resolver import Resolver
from taxon_names_resolver import taxTree

# EXAMPLE NAMES
terms = ['Homo sapiens', 'Gorilla gorilla', 'Pongo pongo', 'Macca mulatta',
'Mus musculus', 'Ailuropoda melanoleuca', 'Ailurus fulgens',
'Chlorotalpa tytonis', 'Arabidopsis thaliana', 'Bacillus subtilus']

# RESOLVE
resolver = Resolver(terms=terms, datasource="NCBI")
resolver.main()

# CREATE TREE
ranks = resolver.retrieve('classification_path_ranks')
idents = resolver.retrieve('query_name')
lineages = resolver.retrieve('classification_path')
treestring = taxTree(idents, ranks, lineages)

# SAVE
with open('example.tre', 'w') as file:
    file.write(treestring)
```

![tree](https://raw.githubusercontent.com/DomBennett/TaxonNamesResolver/master/example.jpg "An example taxonomic tree")

The tree can be visualised with [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) or via one of these websites [Trex-Online](http://www.trex.uqam.ca/index.php?action=newick) or [ETE](http://etetoolkit.org/treeview/).

Alternatively it can be viewed in python using the following code:
```{python}
# IMPORT
from Bio import Phylo
from cStringIO import StringIO

# PROCESS
tree = Phylo.read(StringIO(treestring), "newick")  # Pass to Phylo
Phylo.draw_ascii(tree)  # Draw as ASCII
```

Or `R` if it has been saved as a .tre:
```{R}
# LIBRARY
library(ape)

# PROCESS
tree <- read.tree('example.tre')
plot(tree)
```

`taxTree` by default uses all the taxonomic levels available in NCBI's taxonomy. To change this use the argument `taxonomy` to select with which ranks you would like to build a tree, e.g.

```{python}
taxonomy = ['family', 'superfamily', 'order', 'class', 'kingdom', 'superkingdom']
treestring = taxTree(idents, ranks, lineages, taxonomy=taxonomy)
```

## Acknowledgements
* Thanks to Lawrence Hudson for coding inspiration
* Thanks to the GNR API and the [Global Names](http://www.globalnames.org/) family for creating this wonderful resource

## Author
Dominic John Bennett (dominic.john.bennett@gmail.com)
