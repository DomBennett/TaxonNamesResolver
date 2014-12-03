#! /usr/bin/env python
"""
Tests for misc tools
"""

import unittest
from taxon_names_resolver import misc_tools as mt

# DATA
idents = ['Homo sapiens', 'Pongo pongo', 'Mus musculus', 'Bacillus subtilus',
          'Gorilla gorilla', 'Ailuropoda melanoleuca', 'Ailurus fulgens',
          'Arabidopsis thaliana', 'Macca mulatta', 'Chlorotalpa tytonis']
ranks = [['superkingdom', '', 'kingdom', '', '', '', '', 'phylum',
          'subphylum', '', 'superclass', '', '', '', '', '', 'class',
          '', '', 'superorder', 'order', 'suborder', 'infraorder',
          'parvorder', 'superfamily', 'family', 'subfamily', 'genus',
          'species'],
         ['superkingdom', '', 'kingdom', '', '', '', '', 'phylum',
          'subphylum', '', 'superclass', '', '', '', '', '', 'class',
          '', '', 'superorder', 'order', 'suborder', 'infraorder',
          'parvorder', 'superfamily', 'family', 'subfamily', 'genus'],
         ['superkingdom', '', 'kingdom', '', '', '', '', 'phylum',
          'subphylum', '', 'superclass', '', '', '', '', '', 'class',
          '', '', 'superorder', '', 'order', 'suborder', '', 'family',
          'subfamily', 'genus', 'subgenus', 'species'],
         ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus',
          'species group', 'species'],
         ['superkingdom', '', 'kingdom',
          '', '', '', '', 'phylum', 'subphylum', '', 'superclass', '',
          '', '', '', '', 'class', '', '', 'superorder', 'order',
          'suborder', 'infraorder', 'parvorder', 'superfamily', 'family',
          'subfamily', 'genus', 'species'],
         ['superkingdom', '', 'kingdom', '', '', '', '', 'phylum',
          'subphylum', '', 'superclass', '', '', '', '', '', 'class',
          '', '', 'superorder', 'order', 'suborder', 'family', 'genus',
          'species'],
         ['superkingdom', '', 'kingdom', '', '', '', '', 'phylum',
          'subphylum', '', 'superclass', '', '', '', '', '', 'class',
          '', '', 'superorder', 'order', 'suborder', 'family', 'genus',
          'species'],
         ['superkingdom', 'kingdom', 'phylum', '', '', '', '', '', '',
          '', '', 'subclass', '', 'order', 'family', 'tribe', 'genus',
          'species'],
         ['superkingdom', '', 'kingdom', '', '', '', '', 'phylum',
          'subphylum', '', 'superclass', '', '', '', '', '', 'class',
          '', '', 'superorder', 'order', 'suborder', 'infraorder',
          'parvorder', 'superfamily', 'family', 'subfamily', 'genus',
          'species'],
         ['superkingdom', '', 'kingdom', '', '', '', '', 'phylum',
          'subphylum', '', 'superclass', '', '', '', '', '', 'class',
          '', '', 'superorder', 'family', 'genus']]
lineages = [['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini',
             'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae',
             'Homininae', 'Homo', 'Homo sapiens'],
            ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini',
             'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae',
             'Ponginae', 'Pongo'],
            ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Euarchontoglires', 'Glires', 'Rodentia',
             'Sciurognathi', 'Muroidea', 'Muridae', 'Murinae', 'Mus', 'Mus',
             'Mus musculus'],
            ['Bacteria', 'Firmicutes', 'Bacilli', 'Bacillales', 'Bacillaceae',
             'Bacillus', 'Bacillus subtilis', 'Bacillus subtilis'],
            ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini',
             'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae',
             'Homininae', 'Gorilla', 'Gorilla gorilla'],
            ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Laurasiatheria', 'Carnivora', 'Caniformia',
             'Ursidae', 'Ailuropoda', 'Ailuropoda melanoleuca'],
            ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Laurasiatheria', 'Carnivora', 'Caniformia',
             'Ailuridae', 'Ailurus', 'Ailurus fulgens'],
            ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Streptophytina',
             'Embryophyta', 'Tracheophyta', 'Euphyllophyta', 'Spermatophyta',
             'Magnoliophyta', '', '', '', '', 'Brassicales', 'Brassicaceae',
             'Camelineae', 'Arabidopsis', 'Arabidopsis thaliana'],
            ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini',
             'Simiiformes', 'Catarrhini', 'Cercopithecoidea',
             'Cercopithecidae', 'Cercopithecinae', 'Macaca', 'Macaca mulatta'],
            ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria',
             'Coelomata', 'Deuterostomia', 'Chordata', 'Craniata',
             'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi',
             'Sarcopterygii', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria',
             'Eutheria', 'Afrotheria', 'Chrysochloridae', 'Chlorotalpa']]


class MiscToolsTestSuite(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_taxref(self):
        testref = mt.TaxRef(ident='Genus species', rank='species')
        testref.change(ident='Family genus', rank='genus')
        self.assertEqual(str(testref), 'Family genus -- genus@3')
        self.assertEqual(testref.counter, 1)
        with self.assertRaises(ValueError):
            testref.ident = 48

    def test_taxdict(self):
        taxonomy = ['species', 'family', 'superkingdom']
        taxdict = mt.TaxDict(idents=idents, ranks=ranks, lineages=lineages,
                             taxonomy=taxonomy)
        with self.assertRaises(IndexError):
            # 3 and more will raise an error, only 3 levels in taxonomy
            taxslice = taxdict._slice(3)
        taxslice = taxdict._slice(2)
        # should be list of tuples of len idents
        self.assertEqual(len(taxslice), len(idents))
        # take taxslice at the superkingdom and group
        # should be two groups: eukaryotes and bacteria
        rung = taxdict._group(taxslice)
        self.assertEqual(len(rung), 2)
        # build a hierarchy for the three leveled taxonomy
        hierarchy = taxdict.hierarchy()
        self.assertTrue(all([e in taxonomy for e in hierarchy.keys()]))

    def test_stringclade(self):
        thing1 = mt.TaxRef(ident='thing1', rank='species')
        thing2 = mt.TaxRef(ident='thing2', rank='genus')
        cladestring = mt.stringClade(taxrefs=[thing1, thing2], name='things',
                                     at=4)
        self.assertEqual(cladestring, '(thing1:3.0,thing2:1.0)things')

    def test_taxtree(self):
        treestring = mt.taxTree(idents=idents, ranks=ranks, lineages=lineages)
        self.assertEqual(treestring, '((((Ailurus_fulgens:9.0,Ailuropoda_\
melanoleuca:9.0)Caniformia_suborder:6.0,((((Gorilla_gorilla:4.0,\
Homo_sapiens:4.0)Homininae_subfamily:1.0,Pongo_pongo:3.0)Hominidae_family\
:2.0,Macca_mulatta:7.0)Catarrhini_parvorder:4.0,Mus_musculus:11.0)\
Euarchontoglires_superorder:4.0,Chlorotalpa_tytonis:13.0)Mammalia_class:\
5.0,Arabidopsis_thaliana:20.0)Eukaryota_superkingdom:1.0,Bacillus_subtilus:\
21.0)life;')

if __name__ == '__main__':
    unittest.main()
