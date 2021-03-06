# -----------------------------------------------------------------------------
# Copyright (c) 2015, The WGS-HGT Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkstemp, mkdtemp
from os import close
from os.path import join

from skbio.util import remove_files
from skbio import TreeNode, TabularMSA, Protein

from benchmark.reformat_input import (join_trees,
                                      trim_gene_tree_leaves,
                                      species_gene_mapping,
                                      remove_branch_lengths,
                                      id_mapper,
                                      reformat_rangerdtl,
                                      reformat_trex,
                                      reformat_riatahgt,
                                      reformat_jane4,
                                      reformat_treepuzzle)


class workflowTests(TestCase):
    """ Test WGS-HGT input reformatting functions """

    def setUp(self):
        """
        """
        # test output can be written to this directory
        self.working_dir = mkdtemp()

        # species tree
        f, self.species_tree_fp = mkstemp(prefix='species_',
                                          suffix='.nwk')
        close(f)
        with open(self.species_tree_fp, 'w') as t:
            t.write(species_tree)

        # species tree 2
        f, self.species_tree_2_fp = mkstemp(prefix='species_2_',
                                            suffix='.nwk')
        close(f)
        with open(self.species_tree_2_fp, 'w') as t:
            t.write(species_tree_2)

        # gene tree 1
        f, self.gene_tree_1_fp = mkstemp(prefix='gene_tree_1_',
                                         suffix='.nwk')
        close(f)
        with open(self.gene_tree_1_fp, 'w') as t:
            t.write(gene_tree_1)

        # gene tree 2
        f, self.gene_tree_2_fp = mkstemp(prefix='gene_tree_2_',
                                         suffix='.nwk')
        close(f)
        with open(self.gene_tree_2_fp, 'w') as t:
            t.write(gene_tree_2)

        # gene tree 3
        f, self.gene_tree_3_fp = mkstemp(prefix='gene_tree_3_',
                                         suffix='.nwk')
        close(f)
        with open(self.gene_tree_3_fp, 'w') as t:
            t.write(gene_tree_3)

        # MSA FASTA 3
        f, self.msa_fa_3_fp = mkstemp(prefix='msa_3_',
                                      suffix='.fa')
        close(f)
        with open(self.msa_fa_3_fp, 'w') as t:
            t.write(msa_fa_3)

        self.files_to_remove = [self.species_tree_fp,
                                self.species_tree_2_fp,
                                self.gene_tree_1_fp,
                                self.gene_tree_2_fp,
                                self.gene_tree_3_fp,
                                self.msa_fa_3_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.working_dir)

    def test_join_trees(self):
        """ Test concatenate Newick trees into one file (species, gene)
        """
        self.output_file = join(self.working_dir, 'output_file.nwk')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        join_trees(gene_tree_1, species_tree, self.output_file)
        with open(self.output_file, 'r') as out_f:
            species_gene_tree_1_obs = out_f.read()
        self.assertEqual(species_gene_tree_1_obs, species_gene_tree_1_exp)

    def test_trim_gene_tree_leaves(self):
        """ Test remove '_GENENAME' from tree leaf names (if exists)
        """
        leaves_exp = ["SE001", "SE002", "SE003", "SE004", "SE005", "SE006",
                      "SE007", "SE008", "SE009", "SE010"]
        gene_tree_2 = TreeNode.read(self.gene_tree_2_fp, format='newick')
        trim_gene_tree_leaves(gene_tree_2)
        leaves_obs = sorted([node.name for node in gene_tree_2.tips()])
        self.assertListEqual(leaves_obs, leaves_exp)

    def test_species_gene_mapping(self):
        """ Test finding the association between species and gene tree leaves
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        mapping_exp = {"SE001": ["SE001_01623", "SE001_04123"],
                       "SE002": ["SE002_01623"], "SE003": ["SE003_01623"],
                       "SE004": ["SE004_01623"],
                       "SE005": ["SE005_01623", "SE005_04123"],
                       "SE006": ["SE006_01623", "SE006_04123"],
                       "SE007": ["SE007_01623"],
                       "SE008": ["SE008_01623", "SE008_04123"],
                       "SE009": ["SE009_01623", "SE009_04123"],
                       "SE010": ["SE010_01623", "SE010_04123"]}
        mapping_obs = species_gene_mapping(gene_tree_1, species_tree)
        self.assertDictEqual(dict(mapping_obs), mapping_exp)

    def test_species_gene_mapping_check_species_labels(self):
        """ Test ValueError raised
        """
        species_tree = TreeNode.read(self.species_tree_2_fp, format='newick')
        gene_tree_3 = TreeNode.read(self.gene_tree_3_fp, format='newick')
        self.assertRaises(ValueError,
                          species_gene_mapping,
                          gene_tree=gene_tree_3,
                          species_tree=species_tree)

    def test_remove_branch_lengths(self):
        """ Test removing branch lengths from tree
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        remove_branch_lengths(tree=species_tree)
        species_tree_exp = (
            "(((((((SE001,SE010),SE008),(SE006,SE009)),SE005),"
            "SE004),SE003),(SE002,SE007));")
        self.assertEqual(species_tree_exp, str(species_tree)[:-1])

    def test_id_mapper(self):
        """ Test id_mapper() functionality
        """
        index = [u'SE001/02297', u'SE002/02297', u'SE003/02297',
                 u'SE004/02297', u'SE005/02297', u'SE006/02297',
                 u'SE008/02297', u'SE009/02297', u'SE010/02297']
        mapping = id_mapper(index)
        mapping_exp = {u'SE008/02297': u'SE008',
                       u'SE003/02297': u'SE003',
                       u'SE009/02297': u'SE009',
                       u'SE005/02297': u'SE005',
                       u'SE001/02297': u'SE001',
                       u'SE006/02297': u'SE006',
                       u'SE004/02297': u'SE004',
                       u'SE010/02297': u'SE010',
                       u'SE002/02297': u'SE002'}
        self.assertDictEqual(mapping_exp, mapping)

    def test_reformat_rangerdtl(self):
        """ Test functionality of reformat_rangerdtl()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nwk")
        reformat_rangerdtl(gene_tree_1,
                           species_tree,
                           output_tree_fp)
        joined_tree_exp = [
            "(((((((SE001,SE010),SE008),(SE006,SE009)),SE005),"
            "SE004),SE003),(SE002,SE007));\n",
            "(((((((SE001_01623,SE010_01623),SE008_01623),(SE006_01623,"
            "SE009_01623)),SE005_01623),SE004_01623),SE003_01623),"
            "((SE002_01623,SE007_01623),((((SE001_04123,SE010_04123),"
            "SE008_04123),(SE006_04123,SE009_04123)),SE005_04123)));\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            joined_tree_act = output_tree_f.readlines()
        self.assertListEqual(joined_tree_exp, joined_tree_act)

    def test_reformat_trex(self):
        """ Test functionality of reformat_trex()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nwk")
        reformat_trex(gene_tree_1,
                      species_tree,
                      output_tree_fp)
        reformat_tree_exp = [
            "(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
            "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
            "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
            "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
            "1.6606127,SE007:0.70000178):1.6331374):1.594016;\n",
            "(((((((SE001:2.1494876,SE010:2.1494876):"
            "3.7761166,SE008:5.9256042):0.2102448,(SE006:"
            "5.2329068,SE009:5.2329068):0.9029422):0.2054233,"
            "SE005:6.3412723):0.3714563,SE004:6.7127286):"
            "0.7293362,SE003:7.4420648):0.2444784,((SE002:"
            "6.0534057,SE007:6.0534057):0.4589905,((((SE001:"
            "2.1494876,SE010:2.1494876):3.7761166,SE008:"
            "5.9256042):0.2102448,(SE006:5.2329068,SE009:"
            "5.2329068):0.9029422):0.2054233,SE005:6.3412723):"
            "0.1711239):1.174147):1.594016;\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)

    def test_reformat_riatahgt(self):
        """ Test functionality of reformat_riatahgt()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nex")
        reformat_riatahgt(gene_tree_1,
                          species_tree,
                          output_tree_fp)
        reformat_tree_exp = [
            "#NEXUS\n", "BEGIN TREES;\n",
            "Tree speciesTree = "
            "(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
            "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
            "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
            "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
            "1.6606127,SE007:0.70000178):1.6331374):1.594016;\n",
            "Tree geneTree = "
            "(((((((SE001:2.1494876,SE010:2.1494876):"
            "3.7761166,SE008:5.9256042):0.2102448,(SE006:"
            "5.2329068,SE009:5.2329068):0.9029422):0.2054233,"
            "SE005:6.3412723):0.3714563,SE004:6.7127286):"
            "0.7293362,SE003:7.4420648):0.2444784,((SE002:"
            "6.0534057,SE007:6.0534057):0.4589905,((((SE001:"
            "2.1494876,SE010:2.1494876):3.7761166,SE008:"
            "5.9256042):0.2102448,(SE006:5.2329068,SE009:"
            "5.2329068):0.9029422):0.2054233,SE005:6.3412723):"
            "0.1711239):1.174147):1.594016;\n",
            "END;\n",
            "BEGIN PHYLONET;\n",
            "RIATAHGT speciesTree {geneTree};\n",
            "END;\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)

    def test_reformat_jane4(self):
        """ Test functionality of reformat_jane4()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_1 = TreeNode.read(self.gene_tree_1_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nex")
        reformat_jane4(gene_tree_1,
                       species_tree,
                       output_tree_fp)
        reformat_tree_exp = [
            "#NEXUS\n", "begin host;\n",
            "tree host = "
            "(((((((SE001,SE010),SE008),(SE006,SE009)),SE005),SE004),SE003),"
            "(SE002,SE007));\n", "\n",
            "endblock;\n", "begin parasite;\n",
            "tree parasite = "
            "(((((((SE001_01623,SE010_01623),SE008_01623),(SE006_01623,"
            "SE009_01623)),SE005_01623),SE004_01623),SE003_01623),"
            "((SE002_01623,SE007_01623),((((SE001_04123,SE010_04123),"
            "SE008_04123),(SE006_04123,SE009_04123)),SE005_04123)));\n", "\n",
            "endblock;\n",
            "begin distribution;\n",
            "Range SE010_01623:SE010, SE010_04123:SE010, SE009_01623:SE009, "
            "SE009_04123:SE009, SE008_01623:SE008, SE008_04123:SE008, "
            "SE007_01623:SE007, SE006_01623:SE006, SE006_04123:SE006, "
            "SE005_01623:SE005, SE005_04123:SE005, SE004_01623:SE004, "
            "SE003_01623:SE003, SE002_01623:SE002, SE001_01623:SE001, "
            "SE001_04123:SE001;\n",
            "endblock;\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)

    def test_reformat_treepuzzle(self):
        """ Test functionality of reformat_treepuzzle()
        """
        species_tree = TreeNode.read(self.species_tree_fp, format='newick')
        gene_tree_3 = TreeNode.read(self.gene_tree_3_fp, format='newick')
        output_tree_fp = join(self.working_dir, "joined_trees.nwk")
        output_msa_phy_fp = join(self.working_dir, "gene_tree_3.phy")
        reformat_treepuzzle(gene_tree_3,
                            species_tree,
                            self.msa_fa_3_fp,
                            output_tree_fp,
                            output_msa_phy_fp)
        reformat_tree_exp = [
            "(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
            "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
            "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
            "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
            "1.6606127,SE007:0.70000178):1.6331374);\n",
            "(((((((SE001:2.1494876,SE010:2.1494876):"
            "3.7761166,SE008:5.9256042):0.2102448,(SE006:"
            "5.2329068,SE009:5.2329068):0.9029422):0.2054233,"
            "SE005:6.3412723):0.3714563,SE004:6.7127286):"
            "0.7293362,SE003:7.4420648):0.2444784,SE002:"
            "7.6865432);\n"]
        with open(output_tree_fp, 'r') as output_tree_f:
            reformat_tree_act = output_tree_f.readlines()
        self.assertListEqual(reformat_tree_exp, reformat_tree_act)
        msa_fa = TabularMSA.read(output_msa_phy_fp, constructor=Protein)
        labels_exp = [u'SE001', u'SE002', u'SE003', u'SE004', u'SE005',
                      u'SE006', u'SE008', u'SE009', u'SE010']
        labels_act = list(msa_fa.index)
        self.assertListEqual(labels_exp, labels_act)


# 10 species
species_tree = ("(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
                "0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):"
                "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
                "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
                "1.6606127,SE007:0.70000178):1.6331374):1.594016;")
# 10 species, 1 with same label
species_tree_2 = ("(((((((SE001:2.1494877,SE010:1.08661):3.7761166,SE008:"
                  "0.86305436):0.21024487,(SE001:0.56704221,SE009:0.5014676):"
                  "0.90294223):0.20542323,SE005:3.0992506):0.37145632,SE004:"
                  "1.8129133):0.72933621,SE003:1.737411):0.24447835,(SE002:"
                  "1.6606127,SE007:0.70000178):1.6331374):1.594016;")
# 10 species, 16 genes (gain)
gene_tree_1 = ("(((((((SE001_01623:2.1494876,SE010_01623:2.1494876):"
               "3.7761166,SE008_01623:5.9256042):0.2102448,(SE006_01623:"
               "5.2329068,SE009_01623:5.2329068):0.9029422):0.2054233,"
               "SE005_01623:6.3412723):0.3714563,SE004_01623:6.7127286):"
               "0.7293362,SE003_01623:7.4420648):0.2444784,((SE002_01623:"
               "6.0534057,SE007_01623:6.0534057):0.4589905,((((SE001_04123:"
               "2.1494876,SE010_04123:2.1494876):3.7761166,SE008_04123:"
               "5.9256042):0.2102448,(SE006_04123:5.2329068,SE009_04123:"
               "5.2329068):0.9029422):0.2054233,SE005_04123:6.3412723):"
               "0.1711239):1.174147):1.594016;")
# 10 species, 10 genes
gene_tree_2 = ("(((((((SE001_00009:2.1494876,SE010_00009:2.1494876):"
               "3.7761166,SE008_00009:5.9256042):0.2102448,(SE006_00009:"
               "5.2329068,SE009_00009:5.2329068):0.9029422):0.2054233,"
               "SE005_00009:6.3412723):0.3714563,SE004_00009:6.7127286):"
               "0.7293362,SE003_00009:7.4420648):0.2444784,(SE002_00009:"
               "6.0534057,SE007_00009:6.0534057):1.6331375):1.594016;")
# 10 species, 9 genes (loss)
gene_tree_3 = ("(((((((SE001_02297:2.1494876,SE010_02297:2.1494876):"
               "3.7761166,SE008_02297:5.9256042):0.2102448,(SE006_02297:"
               "5.2329068,SE009_02297:5.2329068):0.9029422):0.2054233,"
               "SE005_02297:6.3412723):0.3714563,SE004_02297:6.7127286):"
               "0.7293362,SE003_02297:7.4420648):0.2444784,SE002_02297:"
               "7.6865432):1.594016;")
# MSA, 9 genes
msa_fa_3 = """>SE001/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTLRKDI\
SVGIDPVKAKKAANNRNSFSAIYKEWYEHKKQVWSVGYASELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMEMANKARRRCGEVFSYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADAIPAFNKALRTFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFPTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE002/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTQLSSIPKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEWPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWVADWLDEKLE
>SE003/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRMIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEWPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EAMQWWADWLDEKVE
>SE004/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSKITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE005/02297
MLTVKQIEKAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKFPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKNRRRCGEVFRYAIVTGAAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKGLATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMRWWADWLDEKVE
>SE006/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRL\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NGKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE008/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDNKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDNKVE
>SE009/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPMMTLQEARDKAWTARKDI\
SVGIDPVKAKKASSNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRL\
EDRGAMERANKARRRCGEVFRYAIVTGRAKYNPAPDLADAMKGYRKKNFPFLPADQIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVNEFVFAGR\
NGKKKPICENAVLLVIKQIGYEGLESGHGFRHEFSTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
>SE010/02297
MLTVKQIEAAKPKERPYRLLDGNGLYLYVPVSGKKVWQLRYKIDGKEKILTVGKYPLMTLQEARDKAWTLRKDI\
SVGIDPVKAKKASNNNNSFSAIYKEWYEHKKQVWSVGYATELAKMFDDDILPIIGGLEIQDIQPMQLLEVIRRF\
EDRGAMEMANKARRRCGEVFSYAIVTGRAKYNPAPDLADAMKGYRGKNFPFLPADAIPAFNKALATFSGSIVSL\
IATKVLRYTALRTKELRSMLWKNVDFENRIITIDASVMKGRKIHVVPMSDQVVELLTTLSSITKPVSEFVFAGR\
NDKKKPICENAVLLVIKQIGYEGLESGHGFRHEFPTIMNEHEYPADAIEVQLAHANGGSVRGIYNHAQYLDKRR\
EMMQWWADWLDEKVE
"""
# concatenated species and gene trees
species_gene_tree_1_exp = """(((((((SE001:2.1494877,SE010:1.08661):3.7761166\
,SE008:0.86305436):0.21024487,(SE006:0.56704221,SE009:0.5014676):0.90294\
223):0.20542323,SE005:3.0992506):0.37145632,SE004:1.8129133):0.72933621,SE00\
3:1.737411):0.24447835,(SE002:1.6606127,SE007:0.70000178):1.6331374):1.59401\
6;
(((((((SE001_01623:2.1494876,SE010_01623:2.1494876):3.7761166,SE008_01623:5.\
9256042):0.2102448,(SE006_01623:5.2329068,SE009_01623:5.2329068):0.90294\
22):0.2054233,SE005_01623:6.3412723):0.3714563,SE004_01623:6.7127286):0.7293\
362,SE003_01623:7.4420648):0.2444784,((SE002_01623:6.0534057,SE007_01623:6.0\
534057):0.4589905,((((SE001_04123:2.1494876,SE010_04123:2.1494876):3.776\
1166,SE008_04123:5.9256042):0.2102448,(SE006_04123:5.2329068,SE009_0\
4123:5.2329068):0.9029422):0.2054233,SE005_04123:6.3412723):0.171123\
9):1.174147):1.594016;
"""


if __name__ == '__main__':
    main()
