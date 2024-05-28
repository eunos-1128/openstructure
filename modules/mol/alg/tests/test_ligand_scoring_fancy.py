import unittest, os, sys
from functools import lru_cache

import numpy as np

import ost
from ost import io, mol, geom
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.ligand_scoring_base import *
    from ost.mol.alg import ligand_scoring_base
    from ost.mol.alg import ligand_scoring_scrmsd
    from ost.mol.alg import ligand_scoring_lddtpli
except ImportError:
    print("Failed to import ligand_scoring.py. Happens when numpy, scipy or "
          "networkx is missing. Ignoring test_ligand_scoring.py tests.")
    sys.exit(0)


def _GetTestfilePath(filename):
    """Get the path to the test file given filename"""
    return os.path.join('testfiles', filename)


@lru_cache(maxsize=None)
def _LoadMMCIF(filename):
    path = _GetTestfilePath(filename)
    ent = io.LoadMMCIF(path)
    return ent


@lru_cache(maxsize=None)
def _LoadPDB(filename):
    path = _GetTestfilePath(filename)
    ent = io.LoadPDB(path)
    return ent


@lru_cache(maxsize=None)
def _LoadEntity(filename):
    path = _GetTestfilePath(filename)
    ent = io.LoadEntity(path)
    return ent


class TestLigandScoringFancy(unittest.TestCase):

    def setUp(self):
        # Silence expected warnings about ignoring of ligands in binding site
        ost.PushVerbosityLevel(ost.LogLevel.Error)

    def tearDown(self):
        ost.PopVerbosityLevel()

    def test_extract_ligands_mmCIF(self):
        """Test that we can extract ligands from mmCIF files.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        sc = LigandScorer(mdl, trg, None, None)

        self.assertEqual(len(sc.target_ligands),  7)
        self.assertEqual(len(sc.model_ligands), 1)
        self.assertEqual(len([r for r in sc.target.residues if r.is_ligand]), 7)
        self.assertEqual(len([r for r in sc.model.residues if r.is_ligand]), 1)

    def test_extract_ligands_PDB(self):
        """Test that we can extract ligands from PDB files containing HET records.
        """
        trg = _LoadPDB("1R8Q.pdb")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        sc = LigandScorer(mdl, trg, None, None)

        self.assertEqual(len(sc.target_ligands),  7)
        self.assertEqual(len(sc.model_ligands), 1)
        self.assertEqual(len([r for r in sc.target.residues if r.is_ligand]), 7)
        self.assertEqual(len([r for r in sc.model.residues if r.is_ligand]), 1)

    def test_init_given_ligands(self):
        """Test that we can instantiate the scorer with ligands contained in
        the target and model entity and given in a list.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        # Pass entity views
        trg_lig = [trg.Select("rname=MG"), trg.Select("rname=G3D")]
        mdl_lig = [mdl.Select("rname=G3D")]
        sc = LigandScorer(mdl, trg, mdl_lig, trg_lig)

        self.assertEqual(len(sc.target_ligands), 4)
        self.assertEqual(len(sc.model_ligands), 1)
        # IsLigand flag should still be set even on not selected ligands
        self.assertEqual(len([r for r in sc.target.residues if r.is_ligand]), 7)
        self.assertEqual(len([r for r in sc.model.residues if r.is_ligand]), 1)

        # Ensure the residues are not copied
        self.assertEqual(len(sc.target.Select("rname=MG").residues), 2)
        self.assertEqual(len(sc.target.Select("rname=G3D").residues), 2)
        self.assertEqual(len(sc.model.Select("rname=G3D").residues), 1)

        # Pass residue handles
        trg_lig = [trg.FindResidue("F", 1), trg.FindResidue("H", 1)]
        mdl_lig = [mdl.FindResidue("L_2", 1)]
        sc = LigandScorer(mdl, trg, mdl_lig, trg_lig)

        self.assertEqual(len(sc.target_ligands), 2)
        self.assertEqual(len(sc.model_ligands), 1)

        # Ensure the residues are not copied
        self.assertEqual(len(sc.target.Select("rname=ZN").residues), 1)
        self.assertEqual(len(sc.target.Select("rname=G3D").residues), 2)
        self.assertEqual(len(sc.model.Select("rname=G3D").residues), 1)

    def test_init_sdf_ligands(self):
        """Test that we can instantiate the scorer with ligands from separate SDF files.

        In order to setup the ligand SDF files, the following code was used:
        for prefix in [os.path.join('testfiles', x) for x in ["1r8q", "P84080_model_02"]]:
            trg = io.LoadMMCIF("%s.cif.gz" % prefix)
            trg_prot = trg.Select("protein=True")
            io.SavePDB(trg_prot, "%s_protein.pdb.gz" % prefix)
            lig_num = 0
            for chain in trg.chains:
                if chain.chain_type == mol.ChainType.CHAINTYPE_NON_POLY:
                    lig_sel = trg.Select("cname=%s" % chain.name)
                    lig_ent = mol.CreateEntityFromView(lig_sel, False)
                    io.SaveEntity(lig_ent, "%s_ligand_%d.sdf" % (prefix, lig_num))
                    lig_num += 1
        """
        mdl = _LoadPDB("P84080_model_02_nolig.pdb")
        mdl_ligs = [_LoadEntity("P84080_model_02_ligand_0.sdf")]
        trg = _LoadPDB("1r8q_protein.pdb.gz")
        trg_ligs = [_LoadEntity("1r8q_ligand_%d.sdf" % i) for i in range(7)]

        # Pass entities
        sc = LigandScorer(mdl, trg, mdl_ligs, trg_ligs)

        self.assertEqual(len(sc.target_ligands), 7)
        self.assertEqual(len(sc.model_ligands), 1)
        # Ensure we set the is_ligand flag
        self.assertEqual(len([r for r in sc.target.residues if r.is_ligand]), 7)
        self.assertEqual(len([r for r in sc.model.residues if r.is_ligand]), 1)

        # Pass residues
        mdl_ligs_res = [mdl_ligs[0].residues[0]]
        trg_ligs_res = [res for ent in trg_ligs for res in ent.residues]

        sc = LigandScorer(mdl, trg, mdl_ligs_res, trg_ligs_res)

        self.assertEqual(len(sc.target_ligands), 7)
        self.assertEqual(len(sc.model_ligands), 1)

    def test_init_reject_duplicate_ligands(self):
        """Test that we reject input if multiple ligands with the same chain
         name/residue number are given.
        """
        mdl = _LoadPDB("P84080_model_02_nolig.pdb")
        mdl_ligs = [_LoadEntity("P84080_model_02_ligand_0.sdf")]
        trg = _LoadPDB("1r8q_protein.pdb.gz")
        trg_ligs = [_LoadEntity("1r8q_ligand_%d.sdf" % i) for i in range(7)]

        # Reject identical model ligands
        with self.assertRaises(RuntimeError):
            sc = LigandScorer(mdl, trg, [mdl_ligs[0], mdl_ligs[0]], trg_ligs)

        # Reject identical target ligands
        lig0 = trg_ligs[0].Copy()
        lig1 = trg_ligs[1].Copy()
        ed1 = lig1.EditXCS()
        ed1.RenameChain(lig1.chains[0], lig0.chains[0].name)
        ed1.SetResidueNumber(lig1.residues[0], lig0.residues[0].number)
        with self.assertRaises(RuntimeError):
            sc = LigandScorer(mdl, trg, mdl_ligs, [lig0, lig1])

    def test__ResidueToGraph(self):
        """Test that _ResidueToGraph works as expected
        """
        mdl_lig = _LoadEntity("P84080_model_02_ligand_0.sdf")

        graph = ligand_scoring_base._ResidueToGraph(mdl_lig.residues[0])
        self.assertEqual(len(graph.edges), 34)
        self.assertEqual(len(graph.nodes), 32)
        # Check an arbitrary node
        self.assertEqual([a for a in graph.adj["14"].keys()], ["13", "29"])

        graph = ligand_scoring_base._ResidueToGraph(mdl_lig.residues[0], by_atom_index=True)
        self.assertEqual(len(graph.edges), 34)
        self.assertEqual(len(graph.nodes), 32)
        # Check an arbitrary node
        self.assertEqual([a for a in graph.adj[13].keys()], [12, 28])

    def test__ComputeSymmetries(self):
        """Test that _ComputeSymmetries works.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        trg_mg1 = trg.FindResidue("E", 1)
        trg_g3d1 = trg.FindResidue("F", 1)
        trg_afb1 = trg.FindResidue("G", 1)
        trg_g3d2 = trg.FindResidue("J", 1)
        mdl_g3d = mdl.FindResidue("L_2", 1)

        sym = ligand_scoring_base.ComputeSymmetries(mdl_g3d, trg_g3d1)
        self.assertEqual(len(sym), 72)

        sym = ligand_scoring_base.ComputeSymmetries(mdl_g3d, trg_g3d1, by_atom_index=True)
        self.assertEqual(len(sym), 72)

        # Test that we can match ions read from SDF
        sdf_lig = _LoadEntity("1r8q_ligand_0.sdf")
        sym = ligand_scoring_base.ComputeSymmetries(trg_mg1, sdf_lig.residues[0], by_atom_index=True)
        self.assertEqual(len(sym), 1)

        # Test that it works with views and only consider atoms in the view
        # Skip PA, PB and O[1-3]A and O[1-3]B in target and model
        # We assume atom index are fixed and won't change
        trg_g3d1_sub_ent = trg_g3d1.Select("aindex>6019")
        trg_g3d1_sub = trg_g3d1_sub_ent.residues[0]
        mdl_g3d_sub_ent = mdl_g3d.Select("aindex>1447")
        mdl_g3d_sub = mdl_g3d_sub_ent.residues[0]

        sym = ligand_scoring_base.ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub)
        self.assertEqual(len(sym), 6)

        sym = ligand_scoring_base.ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub, by_atom_index=True)
        self.assertEqual(len(sym), 6)

        # Substructure matches
        sym = ligand_scoring_base.ComputeSymmetries(mdl_g3d, trg_g3d1_sub, substructure_match=True)
        self.assertEqual(len(sym), 6)

        # Missing atoms only allowed in target, not in model
        with self.assertRaises(NoSymmetryError):
            ligand_scoring_base.ComputeSymmetries(mdl_g3d_sub, trg_g3d1, substructure_match=True)

    def test_SCRMSD(self):
        """Test that SCRMSD works.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        trg_mg1 = trg.FindResidue("E", 1)
        trg_g3d1 = trg.FindResidue("F", 1)
        trg_afb1 = trg.FindResidue("G", 1)
        trg_g3d2 = trg.FindResidue("J", 1)
        mdl_g3d = mdl.FindResidue("L_2", 1)

        rmsd = ligand_scoring_scrmsd.SCRMSD(mdl_g3d, trg_g3d1)
        self.assertAlmostEqual(rmsd, 2.21341e-06, 10)
        rmsd = ligand_scoring_scrmsd.SCRMSD(mdl_g3d, trg_g3d2)
        self.assertAlmostEqual(rmsd, 61.21325, 4)

        # Ensure we raise a NoSymmetryError if the ligand is wrong
        with self.assertRaises(NoSymmetryError):
            ligand_scoring_scrmsd.SCRMSD(mdl_g3d, trg_mg1)
        with self.assertRaises(NoSymmetryError):
            ligand_scoring_scrmsd.SCRMSD(mdl_g3d, trg_afb1)

        # Assert that transform works
        trans = geom.Mat4(-0.999256, 0.00788487, -0.0377333, -15.4397,
                          0.0380652, 0.0473315, -0.998154, 29.9477,
                          -0.00608426, -0.998848, -0.0475963, 28.8251,
                          0, 0, 0, 1)
        rmsd = ligand_scoring_scrmsd.SCRMSD(mdl_g3d, trg_g3d2, transformation=trans)
        self.assertAlmostEqual(rmsd, 0.293972, 5)

        # Assert that substructure matches work
        trg_g3d1_sub = trg_g3d1.Select("aindex>6019").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        # mdl_g3d_sub = mdl_g3d.Select("aindex>1447").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        with self.assertRaises(NoIsomorphicSymmetryError):
            ligand_scoring_scrmsd.SCRMSD(mdl_g3d, trg_g3d1_sub)  # no full match

        # But partial match is OK
        rmsd = ligand_scoring_scrmsd.SCRMSD(mdl_g3d, trg_g3d1_sub, substructure_match=True)
        self.assertAlmostEqual(rmsd, 2.2376232209353475e-06, 8)

        # Ensure it doesn't work the other way around - ie incomplete model is invalid
        with self.assertRaises(NoSymmetryError):
            ligand_scoring_scrmsd.SCRMSD(trg_g3d1_sub, mdl_g3d)  # no full match


    def test_compute_rmsd_scores(self):
        """Test that _compute_scores works.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")
        mdl_lig = io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, [mdl_lig], None)

        # Note: expect warning about Binding site of H.ZN1 not mapped to the model
        self.assertEqual(sc.score_matrix.shape, (7, 1))
        np.testing.assert_almost_equal(sc.score_matrix, np.array(
            [[np.nan],
            [0.04244993],
            [np.nan],
            [np.nan],
            [np.nan],
            [0.29399303],
            [np.nan]]), decimal=5)

    def test_compute_lddtpli_scores(self):
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")
        mdl_lig = io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl_lig], None,
                                                  add_mdl_contacts = False,
                                                  lddt_pli_binding_site_radius = 4.0)
        self.assertEqual(sc.score_matrix.shape, (7, 1))
        self.assertTrue(np.isnan(sc.score_matrix[0, 0]))
        self.assertAlmostEqual(sc.score_matrix[1, 0], 0.99843, 5)
        self.assertTrue(np.isnan(sc.score_matrix[2, 0]))
        self.assertTrue(np.isnan(sc.score_matrix[3, 0]))
        self.assertTrue(np.isnan(sc.score_matrix[4, 0]))
        self.assertAlmostEqual(sc.score_matrix[5, 0], 1.0)
        self.assertTrue(np.isnan(sc.score_matrix[6, 0]))

    def test_check_resnames(self):
        """Test that the check_resname argument works.

        When set to True, it should raise an error if any residue in the
        representation of the binding site in the model has a different
        name than in the reference. Here we manually modify a residue
        name to achieve that effect. This is only relevant for the LDDTPLIScorer
        """
        trg_4c0a = _LoadMMCIF("4c0a.cif.gz")
        trg = trg_4c0a.Select("cname=C or cname=I")

        # Here we modify the name of a residue in 4C0A (THR => TPO in C.15)
        # This residue is in the binding site and should trigger the error
        mdl = ost.mol.CreateEntityFromView(trg, include_exlusive_atoms=False)
        ed = mdl.EditICS()
        ed.RenameResidue(mdl.FindResidue("C", 15), "TPO")
        ed.UpdateXCS()

        with self.assertRaises(RuntimeError):
            sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl.FindResidue("I", 1)], [trg.FindResidue("I", 1)], check_resnames=True)
            sc._compute_scores()

        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl.FindResidue("I", 1)], [trg.FindResidue("I", 1)], check_resnames=False)
        sc._compute_scores()

    def test_added_mdl_contacts(self):

        # binding site for ligand in chain G consists of chains A and B
        prot = _LoadMMCIF("1r8q.cif.gz").Copy()

        # model has the full binding site
        mdl = mol.CreateEntityFromView(prot.Select("cname=A,B,G"), True)

        # chain C has same sequence as chain A but is not in contact
        # with ligand in chain G
        # target has thus incomplete binding site only from chain B
        trg = mol.CreateEntityFromView(prot.Select("cname=B,C,G"), True)

        # if added model contacts are not considered, the incomplete binding
        # site only from chain B is perfectly reproduced by model which also has
        # chain B
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, add_mdl_contacts=False)
        self.assertAlmostEqual(sc.score_matrix[0,0], 1.0, 5)

        # if added model contacts are considered, contributions from chain B are
        # perfectly reproduced but all contacts of ligand towards chain A are
        # added as penalty
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, add_mdl_contacts=True)

        lig = prot.Select("cname=G")
        A_count = 0
        B_count = 0
        for a in lig.atoms:
            close_atoms = mdl.FindWithin(a.GetPos(), sc.lddt_pli_radius)
            for ca in close_atoms:
                cname = ca.GetChain().GetName()
                if cname == "G":
                    pass # its a ligand atom...
                elif cname == "A":
                    A_count += 1
                elif cname == "B":
                    B_count += 1
                
        self.assertAlmostEqual(sc.score_matrix[0,0],
                               B_count/(A_count + B_count), 5)

        # Same as before but additionally we remove residue TRP.66 
        # from chain C in the target to test mapping magic...
        # Chain C is NOT in contact with the ligand but we only
        # add contacts from chain A as penalty that are mappable
        # to the closest chain with same sequence. That would be
        # chain C
        query = "cname=B,G or (cname=C and rnum!=66)"
        trg = mol.CreateEntityFromView(prot.Select(query), True)
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, add_mdl_contacts=True)

        TRP66_count = 0
        for a in lig.atoms:
            close_atoms = mdl.FindWithin(a.GetPos(), sc.lddt_pli_radius)
            for ca in close_atoms:
                cname = ca.GetChain().GetName()
                if cname == "A" and ca.GetResidue().GetNumber().GetNum() == 66:
                    TRP66_count += 1

        self.assertEqual(TRP66_count, 134)

        # remove TRP66_count from original penalty
        self.assertAlmostEqual(sc.score_matrix[0,0],
                               B_count/(A_count + B_count - TRP66_count), 5)        

        # Move a random atom in the model from chain B towards the ligand center
        # chain B is also present in the target and interacts with the ligand,
        # but that atom would be far away and thus adds to the penalty. Since
        # the ligand is small enough, the number of added contacts should be
        # exactly the number of ligand atoms.
        mdl_ed = mdl.EditXCS()
        at = mdl.FindResidue("B", mol.ResNum(8)).FindAtom("NZ")
        mdl_ed.SetAtomPos(at, lig.geometric_center)
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, add_mdl_contacts=True)

        # compared to the last assertAlmostEqual, we add the number of ligand
        # atoms as additional penalties
        self.assertAlmostEqual(sc.score_matrix[0,0],
                               B_count/(A_count + B_count - TRP66_count + \
                               lig.GetAtomCount()), 5)


if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_ligand_scoring.py tests.')