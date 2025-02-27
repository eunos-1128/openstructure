import unittest, os, sys
import ost
from ost import io, mol, geom, seq
# check if we can import: fails if numpy or scipy not available
try:
  from ost.mol.alg.scoring import *
except ImportError:
  print("Failed to import scoring.py. Happens when numpy or scipy "\
        "missing. Ignoring test_scoring.py tests.")
  sys.exit(0)

def _LoadFile(file_name):
  """Helper to avoid repeating input path over and over."""
  return io.LoadPDB(os.path.join('testfiles', file_name))

class TestScorer(unittest.TestCase):

  # compare to hardcoded values - no in depth testing
  # this should be sufficient to flag issues when changes are introduced

  def test_scorer_lddt(self):

    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")

    scorer = Scorer(mdl, trg)

    # check global lDDT values
    self.assertAlmostEqual(scorer.lddt, 0.539, 3)
    self.assertAlmostEqual(scorer.bb_lddt, 0.622, 3)
    self.assertAlmostEqual(scorer.ilddt, 0.282, 3)

    for ch in scorer.model.chains:
      self.assertTrue(ch.name in scorer.local_lddt)
      self.assertEqual(len(scorer.local_lddt[ch.name]), ch.GetResidueCount())
      self.assertTrue(ch.name in scorer.bb_local_lddt)
      self.assertEqual(len(scorer.bb_local_lddt[ch.name]), ch.GetResidueCount())

    for ch in scorer.model.chains:
      self.assertTrue(ch.name in scorer.aa_local_lddt)
      self.assertEqual(len(scorer.aa_local_lddt[ch.name]), ch.GetResidueCount())

    # check some random per-residue/per-atom scores
    self.assertEqual(scorer.local_lddt["B"][mol.ResNum(42)], 0.659, 3)
    self.assertEqual(scorer.local_lddt["A"][mol.ResNum(142)], 0.849, 3)
    self.assertEqual(scorer.bb_local_lddt["B"][mol.ResNum(42)], 0.782, 3)
    self.assertEqual(scorer.bb_local_lddt["A"][mol.ResNum(142)], 0.910, 3)
    self.assertEqual(scorer.aa_local_lddt["B"][mol.ResNum(42)]["CA"], 0.718, 3)
    self.assertEqual(scorer.aa_local_lddt["A"][mol.ResNum(142)]["CB"], 0.837, 3)

    # test stereochemistry checks related behaviour

    # stereochemistry issue in mdl sidechain
    bad_mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    ed = bad_mdl.EditXCS()
    at = bad_mdl.FindResidue("B", mol.ResNum(42)).FindAtom("CD")
    pos = at.GetPos()
    new_pos = geom.Vec3(pos[0], pos[1], pos[2] + 1.0)
    ed.SetAtomPos(at, new_pos)
    scorer = Scorer(bad_mdl, trg)
    # original score without stereochemistry issues: 0.659
    # now it should be much worse but above zero
    penalized_score = scorer.local_lddt["B"][mol.ResNum(42)]
    self.assertTrue(penalized_score < 0.5 and penalized_score > 0.1)

    # let's make it really bad, i.e. involve a backbone atom
    at = bad_mdl.FindResidue("B", mol.ResNum(42)).FindAtom("CA")
    pos = at.GetPos()
    new_pos = geom.Vec3(pos[0], pos[1], pos[2] + 1.0)
    ed.SetAtomPos(at, new_pos)
    scorer = Scorer(bad_mdl, trg)
    # original score without stereochemistry issues: 0.659
    # now it should be 0.0
    penalized_score = scorer.local_lddt["B"][mol.ResNum(42)]
    self.assertEqual(penalized_score, 0.0)

    # let's fiddle around in trg stereochemistry
    bad_trg = _LoadFile("1eud_ref.pdb")
    ed = bad_trg.EditXCS()
    # thats the atom that should map to B.42 in mdl
    at = bad_trg.FindResidue("B", mol.ResNum(5)).FindAtom("CD")
    pos = at.GetPos()
    new_pos = geom.Vec3(pos[0], pos[1], pos[2] + 1.0)
    ed.SetAtomPos(at, new_pos)
    scorer = Scorer(mdl, bad_trg)
    # there is no reference info anymore on the whole sidechain
    # The scores of the sidechain atoms should thus be None
    self.assertTrue(scorer.aa_local_lddt["B"][mol.ResNum(42)]["CD"] is None)
    # but not the sidechain atom scores from backbone atoms!
    self.assertFalse(scorer.aa_local_lddt["B"][mol.ResNum(42)]["CA"] is None)
    # also the full per-residue score is still a valid number
    self.assertFalse(scorer.local_lddt["B"][mol.ResNum(42)] is None)

    bad_trg = _LoadFile("1eud_ref.pdb")
    ed = bad_trg.EditXCS()
    at = bad_trg.FindResidue("B", mol.ResNum(5)).FindAtom("CA")
    pos = at.GetPos()
    new_pos = geom.Vec3(pos[0], pos[1], pos[2] + 1.0)

    ed.SetAtomPos(at, new_pos)
    scorer = Scorer(mdl, bad_trg)

    # all scores should be None now for this residue
    self.assertTrue(scorer.aa_local_lddt["B"][mol.ResNum(42)]["CD"] is None)
    self.assertTrue(scorer.aa_local_lddt["B"][mol.ResNum(42)]["CA"] is None)
    self.assertTrue(scorer.local_lddt["B"][mol.ResNum(42)] is None)


  def test_scorer_qsscore(self):

    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")

    scorer = Scorer(mdl, trg)

    # check qs-score related values
    self.assertAlmostEqual(scorer.qs_global, 0.321, 3)
    self.assertAlmostEqual(scorer.qs_best, 0.932, 3)
    self.assertEqual(len(scorer.qs_interfaces), 1)
    self.assertEqual(scorer.qs_interfaces[0], ("A", "B", "A", "B"))

    # should be equal global scores since we're only dealing with
    # a single interface
    self.assertAlmostEqual(scorer.per_interface_qs_global[0], 0.321, 3)
    self.assertAlmostEqual(scorer.per_interface_qs_best[0], 0.932, 3)

  def test_scorer_rigid_scores(self):
    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")
    scorer = Scorer(mdl, trg)
    self.assertAlmostEqual(scorer.gdtts, 0.616, 3)
    self.assertAlmostEqual(scorer.gdtha, 0.473, 3)
    self.assertAlmostEqual(scorer.rmsd, 2.944, 3)

  def test_scorer_contacts(self):
    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")
    scorer = Scorer(mdl, trg)

    self.assertEqual(len(scorer.model_contacts), 48)
    self.assertEqual(len(scorer.native_contacts), 140)

    self.assertAlmostEqual(scorer.ics, 0.415, 3)
    self.assertAlmostEqual(scorer.ics_precision, 0.812, 3)
    self.assertAlmostEqual(scorer.ics_recall, 0.279, 3)
    self.assertAlmostEqual(scorer.ips, 0.342, 3)
    self.assertAlmostEqual(scorer.ips_precision, 0.891, 3)
    self.assertAlmostEqual(scorer.ips_recall, 0.357, 3)


    # per interface scores should be equal since we're only dealing with one
    # interface
    self.assertEqual(len(scorer.per_interface_ics), 1)
    self.assertAlmostEqual(scorer.per_interface_ics[0], 0.415, 3)
    self.assertAlmostEqual(scorer.per_interface_ics_precision[0], 0.812, 3)
    self.assertAlmostEqual(scorer.per_interface_ics_recall[0], 0.279, 3)

    self.assertEqual(len(scorer.per_interface_ips), 1)
    self.assertAlmostEqual(scorer.per_interface_ips[0], 0.342, 3)
    self.assertAlmostEqual(scorer.per_interface_ips_precision[0], 0.891, 3)
    self.assertAlmostEqual(scorer.per_interface_ips_recall[0], 0.357, 3)

  def test_scorer_dockq(self):
    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")
    scorer = Scorer(mdl, trg)
    self.assertEqual(scorer.dockq_interfaces, [("A", "B", "A", "B")])
    self.assertAlmostEqual(scorer.dockq_scores[0], 0.559, 3)
    self.assertAlmostEqual(scorer.fnat[0], 0.279, 3)
    self.assertAlmostEqual(scorer.fnonnat[0], 0.188, 2)
    self.assertAlmostEqual(scorer.irmsd[0], 0.988, 3)
    self.assertAlmostEqual(scorer.lrmsd[0], 5.533, 3)
    self.assertEqual(scorer.nnat[0], 140)
    self.assertEqual(scorer.nmdl[0], 48)

    # ave and wave values are the same as scorer.dockq_scores[0]
    self.assertAlmostEqual(scorer.dockq_wave, scorer.dockq_scores[0], 7)
    self.assertAlmostEqual(scorer.dockq_ave, scorer.dockq_scores[0], 7)

  def test_scorer_patch_scores(self):
    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")
    scorer = Scorer(mdl, trg)

    patch_qs = scorer.patch_qs
    patch_dockq = scorer.patch_dockq

    # check some random values
    self.assertAlmostEqual(patch_qs["B"][5], 0.925, 3)
    self.assertAlmostEqual(patch_qs["B"][17], 0.966, 3)
    self.assertAlmostEqual(patch_qs["A"][4], 0.973, 3)

    self.assertAlmostEqual(patch_dockq["B"][5], 0.858, 3)
    self.assertAlmostEqual(patch_dockq["B"][17], 0.965, 3)
    self.assertAlmostEqual(patch_dockq["A"][4], 0.973, 3)

  def test_scorer_trimmed_contacts(self):
    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_mdl_partial-dimer.pdb")

    scorer = Scorer(mdl, trg)

    # mdl and trg are the same
    self.assertEqual(scorer.ics, 1.0)
    self.assertEqual(scorer.ips, 1.0)

    # let's cut a critical interface loop in the target
    trg = trg.Select("cname=A or (cname=B and rnum!=157:168)")
    scorer = Scorer(mdl, trg)

    # mdl is now longer than trg which lowers ICS/IPS
    self.assertTrue(scorer.ics < 0.75)
    self.assertTrue(scorer.ips < 0.75)

    # but if we use the trimmed versions, it should go up to 1.0
    # again
    self.assertEqual(scorer.ics_trimmed, 1.0)
    self.assertEqual(scorer.ips_trimmed, 1.0)

    # lets see if the trimmed model has the right
    # residues missing
    for r in scorer.model.residues:
      cname = r.GetChain().GetName()
      rnum = r.GetNumber()
      trimmed_r = scorer.trimmed_model.FindResidue(cname, rnum)
      if cname == "B" and (rnum.num >= 157 and rnum.num <= 168):
        self.assertFalse(trimmed_r.IsValid())
      else:
        self.assertTrue(trimmed_r.IsValid())

  def test_scorer_tmscore(self):
    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")
    scorer = Scorer(mdl, trg)
    self.assertAlmostEqual(scorer.tm_score, 0.711, 3)

  def test_scorer_seqres(self):
    mdl = _LoadFile("1eud_mdl_partial-dimer.pdb")
    trg = _LoadFile("1eud_ref.pdb")
    scorer = Scorer(mdl, trg)
    self.assertTrue(scorer.chain_mapper.seqres is None)
    self.assertTrue(scorer.chain_mapper.trg_seqres_mapping is None)

    seqres = seq.CreateSequenceList()
    seqres.AddSequence(seq.CreateSequence("1", "MAIRHCSYTASRKHLYVDKNTKVICQGFTGK"
                                               "QGTFHSQQALEYGTNLVGGTTPGKGGKTHLGL"
                                               "PVFNTVKEAKEQTGATASVIYVPPPFAAAAIN"
                                               "EAIDAEVPLVVCITEGIPQQDMVRVKHRLLRQ"
                                               "GKTRLIGPNCPGVINPGECKIGIMPGHIHKKG"
                                               "RIGIVSRSGTLTYEAVHQTTQVGLGQSLCVGI"
                                               "GGDPFNGTDFTDCLEIFLNDPATEGIILIGEI"
                                               "GGNAEENAAEFLKQHNSGPKSKPVVSFIAGLT"
                                               "APPGRRMGHAGAIIAGGKGGAKEKITALQSAG"
                                               "VVVSMSPAQLGTTIYKEFEKRKML"))
    seqres.AddSequence(seq.CreateSequence("2", "MVNLQEYQSKKLMSDNGVKVQRFFVADTANEA"
                                               "LEAAKRLNAKEIVLKAQILAGGRGKGVFSSGL"
                                               "KGGVHLTKDPEVVGQLAKQMIGYNLATKQTPK"
                                               "EGVKVNKVMVAEALDISRETYLAILMDRSCNG"
                                               "PVLVGSPQGGVDIEEVAASNPELIFKEQIDII"
                                               "EGIKDSQAQRMAENLGFLGPLQNQAADQIKKL"
                                               "YNLFLKIDATQVEVNPFGETPEGQVVCFDAKI"
                                               "NFDDNAEFRQKDIFAMDDKSENEPIENEAAKY"
                                               "DLKYIGLDGNIACFVNGAGLAMATCDIIFLNG"
                                               "GKPANFLDLGGGVKESQVYQAFKLLTADPKVE"
                                               "AILVNIFGGIVNCAIIANGITKACRELELKVP"
                                               "LVVRLEGTNVHEAQNILTNSGLPITSAVDLED"
                                               "AAKKAVASVTKK"))

    trg_seqres_mapping = {"A": "1",
                          "B": "2"}

    # we simply check if the parameters are propagated to the chain mapper
    scorer = Scorer(mdl, trg, seqres=seqres,
                    trg_seqres_mapping=trg_seqres_mapping,
                    resnum_alignments=True)
    self.assertFalse(scorer.chain_mapper.seqres is None)
    self.assertFalse(scorer.chain_mapper.trg_seqres_mapping is None)


if __name__ == "__main__":
  from ost import testutils
  if testutils.DefaultCompoundLibIsSet():
    testutils.RunTests()
  else:
    print('No compound lib available. Ignoring test_scoring.py tests.')
