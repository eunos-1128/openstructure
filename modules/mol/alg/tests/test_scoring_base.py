import unittest
import os

import ost
from ost import io
from ost import seq
from ost.mol.alg import scoring_base

def _GetTestfilePath(filename):
    """Get the path to the test file given filename"""
    return os.path.join('testfiles', filename)

class TestLigandScoringFancy(unittest.TestCase):

    def test_MMCIFPrep(self):

        poly_ent = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"))
        self.assertEqual(poly_ent.GetChainCount(), 4)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["A", "B", "C", "D"])


        # test enabling extract_nonpoly flag
        poly_ent, non_poly_entities = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                                             extract_nonpoly=True)
        self.assertEqual(poly_ent.GetChainCount(), 4)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["A", "B", "C", "D"])
        self.assertEqual(len(non_poly_entities), 7)
        nonpoly_names = [ent.residues[0].name for ent in non_poly_entities]
        self.assertEqual(nonpoly_names, ["MG", "G3D", "AFB", "ZN", "MG", "G3D", "AFB"])

        # test enabling extract_seqres_mapping flag
        poly_ent, seqres, trg_seqres_mapping = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                                                      extract_seqres_mapping=True)
        self.assertEqual(poly_ent.GetChainCount(), 4)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["A", "B", "C", "D"])

        self.assertEqual(len(seqres), 2)
        seqres_1 = seq.CreateSequence("1", "MGNIFANLFKGLFGKKEMRILMVGLDAAGKTTIL"
                                           "YKLKLGEIVTTIPTIGFNVETVEYKNISFTVWDV"
                                           "GGQDKIRPLWRHYFQNTQGLIFVVDSNDRERVNE"
                                           "AREELMRMLAEDELRDAVLLVFANKQDLPNAMNA"
                                           "AEITDKLGLHSLRHRNWYIQATCATSGDGLYEGL"
                                           "DWLSNQLRNQK")
        seqres_2 = seq.CreateSequence("2", "LEANEGSKTLQRNRKMAMGRKKFNMDPKKGIQFL"
                                           "VENELLQNTPEEIARFLYKGEGLNKTAIGDYLGE"
                                           "REELNLAVLHAFVDLHEFTDLNLVQALRQFLWSF"
                                           "RLPGEAQKIDRMMEAFAQRYCLCNPGVFQSTDTC"
                                           "YVLSYSVIMLNTDLHNPNVRDKMGLERFVAMNRG"
                                           "INEGGDLPEELLRNLYDSIRNEPFKIPEDDGND")

        self.assertEqual(seqres[0].GetName(), seqres_1.GetName())
        self.assertEqual(seqres[0].GetGaplessString(), seqres_1.GetGaplessString())
        self.assertEqual(seqres[1].GetName(), seqres_2.GetName())
        self.assertEqual(seqres[1].GetGaplessString(), seqres_2.GetGaplessString())
        expected_trg_seqres_mapping = {"A": "1",
                                       "B": "2",
                                       "C": "1",
                                       "D": "2"}
        self.assertEqual(len(expected_trg_seqres_mapping), len(trg_seqres_mapping))
        for k,v in trg_seqres_mapping.items():
            self.assertEqual(expected_trg_seqres_mapping[k], v)

        # specify biounit and enable EVERYTHING
        poly_ent, non_poly_entities, seqres, trg_seqres_mapping =\
        scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                               extract_nonpoly=True,
                               extract_seqres_mapping=True,
                               biounit="1")
        self.assertEqual(poly_ent.GetChainCount(), 2)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["1.A", "1.B"])
        self.assertEqual(len(non_poly_entities), 4)
        nonpoly_names = [ent.residues[0].name for ent in non_poly_entities]
        self.assertEqual(nonpoly_names, ["MG", "G3D", "AFB", "ZN"])
        self.assertEqual(seqres[0].GetName(), seqres_1.GetName())
        self.assertEqual(seqres[0].GetGaplessString(), seqres_1.GetGaplessString())
        self.assertEqual(seqres[1].GetName(), seqres_2.GetName())
        self.assertEqual(seqres[1].GetGaplessString(), seqres_2.GetGaplessString())
        expected_trg_seqres_mapping = {"1.A": "1",
                                       "1.B": "2"}
        self.assertEqual(len(expected_trg_seqres_mapping), len(trg_seqres_mapping))
        for k,v in trg_seqres_mapping.items():
            self.assertEqual(expected_trg_seqres_mapping[k], v)

if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_scoring_base.py tests.')