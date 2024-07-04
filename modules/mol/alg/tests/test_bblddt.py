import unittest, os, sys
import ost
from ost import conop
from ost import io, mol, seq, settings
import time
# check if we can import: fails if numpy or scipy not available
try:
    import numpy as np
    from ost.mol.alg.bb_lddt import *
    from ost.mol.alg.lddt import *
    from ost.mol.alg.chain_mapping import *
except ImportError:
    print("Failed to import bb_lddt.py. Happens when numpy or scipy "\
          "missing. Ignoring qsscore.py tests.")
    sys.exit(0)

def _LoadFile(file_name):
    """Helper to avoid repeating input path over and over."""
    return io.LoadPDB(os.path.join('testfiles', file_name))

class TestBBlDDT(unittest.TestCase):

    def test_bblddtentity(self):
        ent = _LoadFile("3l1p.1.pdb")
        ent = BBlDDTEntity(ent)
        self.assertEqual(len(ent.view.chains), 4)
        self.assertEqual(ent.GetChain("A").GetName(), "A")
        self.assertEqual(ent.GetChain("B").GetName(), "B")
        self.assertEqual(ent.GetChain("C").GetName(), "C")
        self.assertEqual(ent.GetChain("D").GetName(), "D")
        self.assertRaises(Exception, ent.GetChain, "E")
        self.assertEqual(ent.chain_names, ["A", "B", "C", "D"])
        self.assertEqual(ent.GetSequence("A"), "DMKALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTISRFEALQLSLKNMSKLRPLLEKWVEEADNNENLQEISKSVQARKRKRTSIENRVRWSLETMFLKSPKPSLQQITHIANQLGLEKDVVRVWFSNRRQKGKR")
        self.assertEqual(ent.GetSequence("B"), "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTISRFEALQLSLKNMSKLRPLLEKWVEEADNNENLQEISKSQARKRKRTSIENRVRWSLETMFLKSPKPSLQQITHIANQLGLEKDVVRVWFSNRRQKGKRS")
        self.assertEqual(ent.GetSequence("C"), "TCCACATTTGAAAGGCAAATGGA")
        self.assertEqual(ent.GetSequence("D"), "ATCCATTTGCCTTTCAAATGTGG")

        # check for a couple of positions with manually extracted values

        # GLU
        pos = ent.GetPos("B")
        self.assertAlmostEqual(pos[5,0], -0.901, places=3)
        self.assertAlmostEqual(pos[5,1], 28.167, places=3)
        self.assertAlmostEqual(pos[5,2], 13.955, places=3)

        # GLY
        pos = ent.GetPos("A")
        self.assertAlmostEqual(pos[23,0], 17.563, places=3)
        self.assertAlmostEqual(pos[23,1], -4.082, places=3)
        self.assertAlmostEqual(pos[23,2], 29.005, places=3)

        # Cytosine
        pos = ent.GetPos("C")
        self.assertAlmostEqual(pos[4,0], 14.796, places=3)
        self.assertAlmostEqual(pos[4,1], 24.653, places=3)
        self.assertAlmostEqual(pos[4,2], 59.318, places=3)


        # check pairwise dist, chain names are always sorted =>
        # A is rows, C is cols 
        dist_one = ent.PairDist("A", "C")
        dist_two = ent.PairDist("C", "A")
        self.assertTrue(np.array_equal(dist_one, dist_two))
        self.assertEqual(dist_one.shape[0], len(ent.GetSequence("A")))
        self.assertEqual(dist_one.shape[1], len(ent.GetSequence("C")))

        # check some random distance between the Gly and Cytosine that we already 
        # checked above
        self.assertAlmostEqual(dist_one[23,4], 41.86, places=2)

        # all chains interact with each other... but hey, check nevertheless
        self.assertEqual(ent.interacting_chains, [("A", "B"), ("A", "C"),
                                                  ("A", "D"), ("B", "C"),
                                                  ("B", "D"), ("C", "D")])

    def test_bb_lddt_scorer(self):

        target = _LoadFile("3l1p.1.pdb")
        model = _LoadFile("3l1p.1_model.pdb")

        # we need to derive a chain mapping prior to scoring
        mapper = ChainMapper(target)
        res = mapper.GetRMSDMapping(model, strategy="greedy_iterative")

        # lets compare with lddt reference implementation

        reference_lddt_scorer = lDDTScorer(target, bb_only=True)

        # make alignments accessible by mdl seq name
        alns = dict()
        for aln in res.alns.values():
            mdl_seq = aln.GetSequence(1)
            alns[mdl_seq.name] = aln

        # lDDT requires a flat mapping with mdl_ch as key and trg_ch as value
        flat_mapping = res.GetFlatMapping(mdl_as_key=True)
        lddt_chain_mapping = dict()
        for mdl_ch, trg_ch in flat_mapping.items():
            if mdl_ch in alns:
                lddt_chain_mapping[mdl_ch] = trg_ch

        reference_lddt_score = reference_lddt_scorer.lDDT(model,
                                                          chain_mapping = lddt_chain_mapping,
                                                          residue_mapping = alns,
                                                          check_resnames=False)[0]

        bb_lddt_scorer = BBlDDTScorer.FromMappingResult(res)
        bb_lddt_score = bb_lddt_scorer.Score(res.mapping)

        self.assertAlmostEqual(reference_lddt_score, bb_lddt_score, places = 4)

if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_bblddt.py tests.')
