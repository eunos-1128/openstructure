import unittest, os, sys
import ost
from ost import conop
from ost import io, mol, seq, settings
from ost.mol.alg import dockq
import json

def _LoadFile(file_name):
    """Helper to avoid repeating input path over and over."""
    return io.LoadPDB(os.path.join('testfiles', file_name))

class TestDockQ(unittest.TestCase):

    def test_dockq(self):

        mdl = _LoadFile("lddtbs_mdl.pdb").Select("peptide=True")
        trg = _LoadFile("lddtbs_ref_1r8q.1.pdb").Select("peptide=True")

        dockq_result = dockq.DockQ(mdl, trg, "A", "B", "B", "A")

        # compare to results computed with DockQ from https://github.com/bjornwallner/DockQ
        # commit: 7a0a1c49ec70e2db68cb160abbe6aaf2f844a4ec
        self.assertEqual(dockq_result["nnat"], 87)
        self.assertEqual(dockq_result["nmdl"], 109)
        self.assertAlmostEqual(dockq_result["fnat"], 0.828, places=3)
        self.assertAlmostEqual(dockq_result["fnonnat"], 0.339, places=3)
        self.assertAlmostEqual(dockq_result["irmsd"], 2.411, places=3)
        self.assertAlmostEqual(dockq_result["lrmsd"], 6.568, places=3)
        self.assertAlmostEqual(dockq_result["DockQ"], 0.578, places=3)

        # enable "peptide mode" by passing -capri_peptide flag to DockQ executable
        # not that we're dealing with a peptide here, but we can still check
        # the slightly different parameterization...
        dockq_result = dockq.DockQ(mdl, trg, "A", "B", "B", "A",
                                   contact_dist_thresh = 4.0,
                                   interface_dist_thresh = 8.0,
                                   interface_cb = True)

        self.assertEqual(dockq_result["nnat"], 54)
        self.assertEqual(dockq_result["nmdl"], 76)
        self.assertAlmostEqual(dockq_result["fnat"], 0.815, places=3)
        self.assertAlmostEqual(dockq_result["fnonnat"], 0.421, places=3)
        self.assertAlmostEqual(dockq_result["irmsd"], 1.579, places=3)
        # for whatever reason, DockQ produces a slightly different number for
        # lrmsd than above... Let's just slightly reduce accuracy
        self.assertAlmostEqual(dockq_result["lrmsd"], 6.569, places=2)
        self.assertAlmostEqual(dockq_result["DockQ"], 0.638, places=3)

if __name__ == "__main__":
    from ost import testutils
    testutils.RunTests()
