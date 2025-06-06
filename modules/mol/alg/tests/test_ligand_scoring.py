import unittest, os, sys
from functools import lru_cache

import numpy as np

import ost
from ost import io, mol, geom, conop, seq
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.ligand_scoring_base import *
    from ost.mol.alg import scoring_base
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
        # they're extracted in the MMCIFPrep function actually
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly=True)
        mdl, mdl_lig = scoring_base.MMCIFPrep(_GetTestfilePath("P84080_model_02.cif.gz"),
                                              extract_nonpoly=True)


        sc = LigandScorer(mdl, trg, mdl_lig, trg_lig)

        self.assertEqual(len(sc.target_ligands),  7)
        self.assertEqual(len(sc.model_ligands), 1)
        self.assertEqual(len([r for r in sc.target_ligands if r.is_ligand]), 7)
        self.assertEqual(len([r for r in sc.model_ligands if r.is_ligand]), 1)

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
        mdl = scoring_base.PDBPrep(_GetTestfilePath("P84080_model_02_nolig.pdb"))
        mdl_ligs = [_LoadEntity("P84080_model_02_ligand_0.sdf")]
        trg = scoring_base.PDBPrep(_GetTestfilePath("1r8q_protein.pdb.gz"))
        trg_ligs = [_LoadEntity("1r8q_ligand_%d.sdf" % i) for i in range(7)]

        # Pass entities
        sc = LigandScorer(mdl, trg, mdl_ligs, trg_ligs)

        self.assertEqual(len(sc.target_ligands), 7)
        self.assertEqual(len(sc.model_ligands), 1)
        # Ensure we set the is_ligand flag
        self.assertEqual(len([r for r in sc.target_ligands if r.is_ligand]), 7)
        self.assertEqual(len([r for r in sc.model_ligands if r.is_ligand]), 1)

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
        mdl = scoring_base.PDBPrep(_GetTestfilePath("P84080_model_02_nolig.pdb"))
        mdl_ligs = [_LoadEntity("P84080_model_02_ligand_0.sdf")]
        trg = scoring_base.PDBPrep(_GetTestfilePath("1r8q_protein.pdb.gz"))
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

    def test_ComputeSymmetries(self):
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
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly=True)
        mdl = scoring_base.MMCIFPrep(_GetTestfilePath("P84080_model_02.cif.gz"))
        mdl_lig = io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))

        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, [mdl_lig], trg_lig)

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
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly=True)

        mdl = scoring_base.MMCIFPrep(_GetTestfilePath("P84080_model_02.cif.gz"))
        mdl_lig = io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))


        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl_lig], trg_lig,
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

    def test_added_mdl_contacts(self):

        # binding site for ligand in chain G consists of chains A and B
        prot = _LoadMMCIF("1r8q.cif.gz").Copy()

        # model has the full binding site
        mdl = mol.CreateEntityFromView(prot.Select("cname=A,B"), True)
        mdl_lig = mol.CreateEntityFromView(prot.Select("cname=G"), True)

        # chain C has same sequence as chain A but is not in contact
        # with ligand in chain G
        # target has thus incomplete binding site only from chain B
        trg = mol.CreateEntityFromView(prot.Select("cname=B,C"), True)
        trg_lig = mol.CreateEntityFromView(prot.Select("cname=G"), True)

        # if added model contacts are not considered, the incomplete binding
        # site only from chain B is perfectly reproduced by model which also has
        # chain B
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl_lig], [trg_lig],
                                                  add_mdl_contacts=False)
        self.assertAlmostEqual(sc.score_matrix[0,0], 1.0, 5)

        # if added model contacts are considered, contributions from chain B are
        # perfectly reproduced but all contacts of ligand towards chain A are
        # added as penalty
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl_lig], [trg_lig],
                                                  add_mdl_contacts=True)

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
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl_lig], [trg_lig],
                                                  add_mdl_contacts=True)

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
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl_lig], [trg_lig],
                                                  add_mdl_contacts=True)

        # compared to the last assertAlmostEqual, we add the number of ligand
        # atoms as additional penalties
        self.assertAlmostEqual(sc.score_matrix[0,0],
                               B_count/(A_count + B_count - TRP66_count + \
                               lig.GetAtomCount()), 5)

    def test_assignment(self):
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly = True)
        mdl, mdl_lig = scoring_base.MMCIFPrep(_GetTestfilePath("P84080_model_02.cif.gz"),
                                              extract_nonpoly = True)
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, trg_lig)
        self.assertEqual(sc.assignment, [(1, 0)])

        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, mdl_lig, trg_lig)
        self.assertEqual(sc.assignment, [(5, 0)])

    def test_dict_results_rmsd(self):
        """Test that the scores are computed correctly
        """
        # 4C0A has more ligands
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                                     extract_nonpoly = True)

        trg_4c0a, trg_4c0a_lig = scoring_base.MMCIFPrep(_GetTestfilePath("4c0a.cif.gz"),
                                                               extract_nonpoly = True)
        sc = ligand_scoring_scrmsd.SCRMSDScorer(trg, trg_4c0a, trg_lig, trg_4c0a_lig)
        expected_keys = {"J", "F"}
        self.assertFalse(expected_keys.symmetric_difference(sc.score.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.aux.keys()))
        # rmsd
        self.assertAlmostEqual(sc.score["J"][mol.ResNum(1)], 0.8016608357429504, 5)
        self.assertAlmostEqual(sc.score["F"][mol.ResNum(1)], 0.9286373257637024, 5)
        # rmsd_details
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["chain_mapping"], {'F': 'D', 'C': 'C'})
        self.assertEqual(len(sc.aux["J"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.aux["J"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.aux["J"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["chain_mapping"], {'B': 'B', 'G': 'A'})
        self.assertEqual(len(sc.aux["F"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.aux["F"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.aux["F"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')

    def test_dict_results_lddtpli(self):
        """Test that the scores are computed correctly
        """
        # 4C0A has more ligands
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly = True)

        trg_4c0a, trg_4c0a_lig = scoring_base.MMCIFPrep(_GetTestfilePath("4c0a.cif.gz"),
                                                        extract_nonpoly = True)

        sc = ligand_scoring_lddtpli.LDDTPLIScorer(trg, trg_4c0a, trg_lig, trg_4c0a_lig,
                                                  add_mdl_contacts=False,
                                                  lddt_pli_binding_site_radius = 4.0)
        expected_keys = {"J", "F"}
        self.assertFalse(expected_keys.symmetric_difference(sc.score.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.aux.keys()))

        # lddt_pli
        self.assertAlmostEqual(sc.score["J"][mol.ResNum(1)], 0.9127105666156202, 5)
        self.assertAlmostEqual(sc.score["F"][mol.ResNum(1)], 0.915929203539823, 5)
        # lddt_pli_details
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["lddt_pli_n_contacts"], 653)
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["lddt_pli_n_contacts"], 678)
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')

        # lddt_pli with added mdl contacts
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(trg, trg_4c0a, trg_lig, trg_4c0a_lig,
                                                  add_mdl_contacts=True)
        self.assertAlmostEqual(sc.score["J"][mol.ResNum(1)], 0.8988340192043895, 5)
        self.assertAlmostEqual(sc.score["F"][mol.ResNum(1)], 0.9039735099337749, 5)
        # lddt_pli_details
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["lddt_pli_n_contacts"], 729)
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.aux["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["lddt_pli_n_contacts"], 755)
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.aux["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')

    def test_ignore_binding_site(self):
        """Test that we ignore non polymer stuff in the binding site.
         NOTE: we should consider changing this behavior in the future and take
         other ligands, peptides and short oligomers into account for superposition.
         When that's the case this test should be adapter
         """
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1SSP.cif.gz"),
                                              extract_nonpoly=True)

        sc = ligand_scoring_scrmsd.SCRMSDScorer(trg, trg, trg_lig, trg_lig)
        expected_bs_ref_res = ['C.GLY62', 'C.GLN63', 'C.ASP64', 'C.PRO65', 'C.TYR66', 'C.CYS76', 'C.PHE77', 'C.ASN123', 'C.HIS187']
        ost.PushVerbosityLevel(ost.LogLevel.Error)
        self.assertEqual([str(r) for r in sc.aux["D"][1]["bs_ref_res"]], expected_bs_ref_res)
        ost.PopVerbosityLevel()

    def test_substructure_match(self):
        """Test that substructure_match=True works."""
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        trg_g3d1 = trg.FindResidue("F", 1)
        mdl_g3d = mdl.FindResidue("L_2", 1)

        # Skip PA, PB and O[1-3]A and O[1-3]B in target and model
        # ie 8 / 32 atoms => coverage 0.75
        # We assume atom index are fixed and won't change
        trg_g3d1_sub_ent = trg_g3d1.Select("aindex>6019")
        trg_g3d1_sub = trg_g3d1_sub_ent.residues[0]

        # without enabling substructure matches
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl.Select("protein=True"), trg.Select("protein=True"),
                                                model_ligands=[mdl_g3d], target_ligands=[trg_g3d1_sub],
                                                substructure_match=False)
        self.assertEqual(sc.coverage_matrix.shape, (1,1))
        self.assertTrue(np.isnan(sc.coverage_matrix[0,0]))
        self.assertEqual(sc.state_matrix[0,0], 3) # error encoding for that particular issue

        # Substructure matches
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl.Select("protein=True"), trg.Select("protein=True"),
                                                model_ligands=[mdl_g3d], target_ligands=[trg_g3d1_sub],
                                                substructure_match=True)
        self.assertEqual(sc.coverage_matrix.shape, (1,1))
        self.assertEqual(sc.coverage_matrix[0,0], 0.75)
        self.assertEqual(sc.state_matrix[0,0], 0) # no error encoded in state

    def test_6jyf(self):
        """6JYF initially caused issues in the CASP15-CAMEO/LIGATE paper where
         the ligand RET was wrongly assigned to short copies of OLA that float
          around and yielded higher scores.
          Here we test that this is resolved correctly."""
        mdl = scoring_base.PDBPrep(_GetTestfilePath("6jyf_mdl.pdb"))
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("6jyf_trg.cif"),
                                              extract_nonpoly=True)
        mdl_lig = _LoadEntity("6jyf_RET_pred.sdf")
        mdl_lig_full = _LoadEntity("6jyf_RET_pred_complete.sdf")

        # Problem is easily fixed by just prioritizing full coverage
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, [mdl_lig], trg_lig,
                                                substructure_match=True)
        self.assertEqual(len(sc.assignment), 1) # only one mdl ligand => 1 assignment
        trg_lig_idx, mdl_lig_idx = sc.assignment[0]
        self.assertEqual(sc.coverage_matrix[trg_lig_idx, mdl_lig_idx], 1.0)
        self.assertEqual(sc.aux['00001_'][1]["target_ligand"].name, "RET")
        self.assertAlmostEqual(sc.score['00001_'][1], 15.56022, 4)
        self.assertAlmostEqual(sc.coverage_matrix[0,0], 1.)
        self.assertTrue(np.isnan(sc.coverage_matrix[1,0]))
        self.assertTrue(np.isnan(sc.coverage_matrix[2,0]))
        self.assertTrue(np.isnan(sc.coverage_matrix[3,0]))
        self.assertTrue(np.isnan(sc.coverage_matrix[4,0]))
        self.assertTrue(np.isnan(sc.coverage_matrix[5,0]))
        self.assertTrue(np.isnan(sc.coverage_matrix[6,0]))
        self.assertTrue(np.isnan(sc.coverage_matrix[7,0]))
        self.assertAlmostEqual(sc.coverage_matrix[8,0], 0.5)
        self.assertAlmostEqual(sc.coverage_matrix[9,0], 0.3)
        self.assertAlmostEqual(sc.coverage_matrix[10,0], 0.45)
        self.assertTrue(np.isnan(sc.coverage_matrix[11,0]))
        self.assertTrue(np.isnan(sc.coverage_matrix[12,0]))
        self.assertAlmostEqual(sc.coverage_matrix[13,0], 0.55)

        # We need to make sure that it also works if the match is partial.
        # For that we load the complete ligand incl. the O missing in target
        # with a coverage of around 95% only.
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, [mdl_lig_full], trg_lig,
                                                substructure_match=True)
        self.assertEqual(len(sc.assignment), 1) # only one mdl ligand => 1 assignment
        trg_lig_idx, mdl_lig_idx = sc.assignment[0]
        self.assertAlmostEqual(sc.coverage_matrix[trg_lig_idx, mdl_lig_idx],0.95238096)
        self.assertEqual(sc.aux['00001_'][1]["target_ligand"].name, "RET")
        self.assertAlmostEqual(sc.score['00001_'][1], 15.56022, 4)

        # Next, we check that coverage_delta has an effect. With a large
        # delta of 0.5 we will assign to OLA which has a higher RMSD
        # but a coverage of 0.52 only.
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, [mdl_lig_full], trg_lig,
                                                substructure_match=True,
                                                coverage_delta=0.5)
        self.assertEqual(len(sc.assignment), 1) # only one mdl ligand => 1 assignment
        trg_lig_idx, mdl_lig_idx = sc.assignment[0]
        self.assertAlmostEqual(sc.coverage_matrix[trg_lig_idx, mdl_lig_idx],  0.52380955)
        self.assertEqual(sc.aux['00001_'][1]["target_ligand"].name, "OLA")
        self.assertAlmostEqual(sc.score['00001_'][1], 6.13006878, 4)

    def test_skip_too_many_symmetries(self):
        """
        Test behaviour of max_symmetries.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        # Pass entity views
        trg_lig = [trg.Select("cname=F")]
        mdl_lig = [mdl.Select("rname=G3D")]

        # G3D has 72 isomorphic mappings to itself.
        # Limit to 10 to raise
        symmetries = ligand_scoring_base.ComputeSymmetries(mdl_lig[0], trg_lig[0], max_symmetries=100)
        self.assertEqual(len(symmetries), 72)
        with self.assertRaises(TooManySymmetriesError):
            ligand_scoring_base.ComputeSymmetries(mdl_lig[0], trg_lig[0], max_symmetries=10)

        # Check the unassignment
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, trg_lig,
                                                max_symmetries=10)

        self.assertFalse("L_2" in sc.score)
        self.assertEqual(sc.assignment, [])
        self.assertEqual(sc.unassigned_target_ligands, [0])
        self.assertEqual(sc.unassigned_model_ligands, [0])

        trg_report, trg_pair_report = sc.get_target_ligand_state_report(0)
        mdl_report, mdl_pair_report = sc.get_model_ligand_state_report(0)

        # the individual ligands are OK
        self.assertEqual(trg_report["short desc"], "OK")
        self.assertEqual(mdl_report["short desc"], "OK")

        # but there are too many symmetries
        self.assertEqual(len(trg_pair_report), 1)
        self.assertEqual(len(mdl_pair_report), 1)
        self.assertEqual(trg_pair_report[0]["short desc"], "symmetries")
        self.assertEqual(mdl_pair_report[0]["short desc"], "symmetries")

    def test_no_binding_site(self):
        """
        Test the behavior when there's no binding site in proximity of
        the ligand. This test was introduced to identify some subtle issues
        with the ligand assignment that can cause it to enter an infinite
        loop when the data matrices are not filled properly.
        """
        trg = _LoadMMCIF("1r8q.cif.gz").Copy()
        

        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly=True)

        mdl = trg.Copy()
        mdl_lig = [l.Copy() for l in trg_lig]

        trg_g3d = trg_lig[1]
        trg_zn = trg_lig[3]
        self.assertEqual(trg_g3d.chains[0].name, "F")       
        self.assertEqual(trg_zn.chains[0].name, "H")
        self.assertEqual(trg_zn.GetAtomCount(), 1)    

        # Move the zinc out of the reference binding site...
        ed = trg_zn.EditXCS()
        ed.SetAtomPos(trg_zn.atoms[0],
                      trg_zn.atoms[0].pos + geom.Vec3(6, 0, 0))

        # Remove some atoms from G3D to decrease coverage. This messed up
        # the assignment in the past.
        ed = trg_g3d.EditXCS()
        ed.DeleteAtom(trg_g3d.residues[0].FindAtom("O6"))
        ed.UpdateICS()

        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, target_ligands=[trg_zn, trg_g3d],
                                                coverage_delta=0, substructure_match=True)

        self.assertTrue(np.isnan(sc.score_matrix[0, 3]))

        trg_report, trg_pair_report = sc.get_target_ligand_state_report(0)

        exp_lig_report = {'state': 10.0,
                          'short desc': 'target_binding_site',
                          'desc': 'No residues were in proximity of the target ligand.'}

        exp_pair_report = [{'state': 1, 'short desc': 'identity',
                            'desc': 'Ligands could not be matched (by subgraph isomorphism)',
                            'indices': [0, 1, 2, 4, 5, 6]},
                           {'state': 6, 'short desc': 'single_ligand_issue',
                            'desc': 'Cannot compute valid pairwise score as either model or target ligand have non-zero state.',
                            'indices': [3]}]

        # order of report is fix
        self.assertDictEqual(trg_report, exp_lig_report)
        self.assertDictEqual(trg_pair_report[0], exp_pair_report[0])
        self.assertDictEqual(trg_pair_report[1], exp_pair_report[1])


    def test_no_lddt_pli_contact(self):
        """
        Test behaviour where a binding site has no LDDT-PLI contacts.

        We give two copies of the target ligand which have binding site atoms
        within radius=5A but no atoms at 4A. We set lddt_pli_radius=4 so that
        there are no contacts for the LDDT-PLI computation, and LDDT is None.

        We check that:
        - We don't get an LDDT-PLI assignment
        - Both target ligands are unassigned and have the
        - We get an RMSD assignment
        - The second copy of the target and model ligands ensure that the
          disambiguation code (which checks for the best LDDT-PLI when 2 RMSDs
          are identical) works in this case (where there is no LDDT-PLI to
          disambiguate the RMSDs).
        - We get LDDT-PLI = None with RMSD assignment
        """
        trg = scoring_base.PDBPrep(_GetTestfilePath("T1118v1.pdb"))
        trg_lig = _LoadEntity("T1118v1_001.sdf")
        mdl = scoring_base.PDBPrep(_GetTestfilePath("T1118v1LG035_1.pdb"))
        mdl_lig = _LoadEntity("T1118v1LG035_1_1_1.sdf")

        # Ensure it's unassigned in LDDT
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [mdl_lig, mdl_lig],
                                                  [trg_lig, trg_lig],
                                                  lddt_pli_radius=4,
                                                  rename_ligand_chain=True)

        # assignment should be empty
        self.assertEqual(len(sc.assignment), 0)
        self.assertEqual(sc.unassigned_target_ligands, [0, 1])
        self.assertEqual(sc.unassigned_model_ligands, [0, 1])

        expected_report = [{'state': 10, 'short desc': 'no_contact',
                            'desc': 'There were no LDDT contacts between the binding site and the ligand, and LDDT-PLI is undefined.',
                            'indices': [0, 1]}]

        # both target ligands are expected to have the same report => no_contact
        report_1, report_pair_1 = sc.get_target_ligand_state_report(0)
        report_2, report_pair_2 = sc.get_target_ligand_state_report(1)
        self.assertEqual(len(report_pair_1), 1)
        self.assertEqual(len(report_pair_2), 1)
        self.assertDictEqual(report_pair_1[0], expected_report[0])
        self.assertDictEqual(report_pair_2[0], expected_report[0])


        # However RMSD should be assigned
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, [mdl_lig, mdl_lig],
                                                [trg_lig, trg_lig],
                                                bs_radius=5,
                                                rename_ligand_chain=True)

        self.assertEqual(sc.assignment, [(0,0), (1,1)])
        self.assertEqual(sc.unassigned_target_ligands, [])
        self.assertEqual(sc.unassigned_model_ligands, [])

        self.assertAlmostEqual(sc.score['00001_FE'][1], 2.1703, 4)
        self.assertAlmostEqual(sc.score['00001_FE_2'][1], 2.1703, 4)

        expected_report = [{'state': 0, 'short desc': 'OK',
                            'desc': 'OK',
                            'indices': [0, 1]}]

        # both target ligands are expected to have the same report => no_contact
        report_1, report_pair_1 = sc.get_target_ligand_state_report(0)
        report_2, report_pair_2 = sc.get_target_ligand_state_report(1)
        self.assertEqual(len(report_pair_1), 1)
        self.assertEqual(len(report_pair_2), 1)
        self.assertDictEqual(report_pair_1[0], expected_report[0])
        self.assertDictEqual(report_pair_2[0], expected_report[0])

    def test_unassigned_reasons(self):
        """Test reasons for being unassigned."""
        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly=True)
        mdl, mdl_lig = scoring_base.MMCIFPrep(_GetTestfilePath("P84080_model_02.cif.gz"),
                                              extract_nonpoly=True)

        # Add interesting ligands to model and target
        mdl_lig_ent = mol.CreateEntity()
        mdl_lig_ed = mdl_lig_ent.EditXCS()
        trg_lig_ent = mol.CreateEntity()
        trg_lig_ed = trg_lig_ent.EditXCS()

        # Add ZN: representation in the model (chain missing in model)
        new_chain = mdl_lig_ed.InsertChain("L_ZN")
        mdl_lig_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        trg_zn_ent = trg_lig[3]
        self.assertEqual(trg_zn_ent.residues[0].name, "ZN")
        new_res = mdl_lig_ed.AppendResidue(new_chain, "ZN")
        new_atom = mdl_lig_ed.InsertAtom(new_res, "ZN",
                                         trg_zn_ent.atoms[0].GetPos(), "ZN")
        new_res.is_ligand = True

        # Add NA: not in contact with target
        new_chain = trg_lig_ed.InsertChain("L_NA")
        trg_lig_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = trg_lig_ed.AppendResidue(new_chain, "NA")
        new_atom = trg_lig_ed.InsertAtom(new_res, "NA", geom.Vec3(100, 100, 100), "NA")
        new_res.is_ligand = True
        new_chain = mdl_lig_ed.InsertChain("L_NA")
        mdl_lig_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = mdl_lig_ed.AppendResidue(new_chain, "NA")
        new_atom = mdl_lig_ed.InsertAtom(new_res, "NA", geom.Vec3(100, 100, 100), "NA")
        new_res.is_ligand = True

        # Add OXY: no symmetry/ not identical -
        new_chain = mdl_lig_ed.InsertChain("L_OXY")
        mdl_lig_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = mdl_lig_ed.AppendResidue(new_chain, "OXY")
        new_atom1 = mdl_lig_ed.InsertAtom(new_res, "O1", geom.Vec3(0, 0, 0), "O")
        new_atom2 = mdl_lig_ed.InsertAtom(new_res, "O2", geom.Vec3(1, 1, 1), "O")
        mdl_lig_ed.Connect(new_atom1, new_atom2)
        new_res.is_ligand = True

        # Add CMO: disconnected
        new_chain = mdl_lig_ed.InsertChain("L_CMO")
        mdl_lig_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = mdl_lig_ed.AppendResidue(new_chain, "CMO")
        new_atom1 = mdl_lig_ed.InsertAtom(new_res, "O", geom.Vec3(0, 0, 0), "O")
        new_atom2 = mdl_lig_ed.InsertAtom(new_res, "C", geom.Vec3(1, 1, 1), "O")
        new_res.is_ligand = True
        new_chain = trg_lig_ed.InsertChain("L_CMO")
        trg_lig_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = trg_lig_ed.AppendResidue(new_chain, "CMO")
        new_atom1 = trg_lig_ed.InsertAtom(new_res, "O", geom.Vec3(0, 0, 0), "O")
        new_atom2 = trg_lig_ed.InsertAtom(new_res, "C", geom.Vec3(1, 1, 1), "O")
        new_res.is_ligand = True

        # Add 3 MG in model: assignment/stoichiometry
        mg_pos = [
            geom.Vec3(3.871, 12.343, 44.485),
            geom.Vec3(3.871, 12.343, 44.485) + 1,
            geom.Vec3(3.871, 12.343, 44.485) + 100
        ]
        for i in range(3):
            new_chain = mdl_lig_ed.InsertChain("L_MG_%d" % i)
            mdl_lig_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
            new_res = mdl_lig_ed.AppendResidue(new_chain, "MG")
            new_atom = mdl_lig_ed.InsertAtom(new_res, "MG", mg_pos[i], "MG")
            new_res.is_ligand = True

        mdl_lig_ed.UpdateICS()
        trg_lig_ed.UpdateICS()

        trg_lig.append(trg_lig_ent)
        mdl_lig.append(mdl_lig_ent)

        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, mdl_lig, trg_lig)

        # Check unassigned targets
        # NA: not in contact with target
        trg_na = sc.target.FindResidue("L_NA", 1)
        self.assertEqual(sc.unassigned_target_ligands_reasons["L_NA"][1], "no_contact")
        # AFB: not identical to anything in the model
        trg_afb = sc.target.FindResidue("G", 1)
        self.assertEqual(sc.unassigned_target_ligands_reasons["G"][1], "identity")
        # F.G3D1: J.G3D1 assigned instead
        trg_fg3d1 = sc.target.FindResidue("F", 1)
        self.assertEqual(sc.unassigned_target_ligands_reasons["F"][1], "stoichiometry")
        # CMO: disconnected
        trg_cmo1 = sc.target.FindResidue("L_CMO", 1)
        self.assertEqual(sc.unassigned_target_ligands_reasons["L_CMO"][1], "disconnected")
        # J.G3D1: assigned to L_2.G3D1 => check if it is assigned
        self.assertTrue(5 not in sc.unassigned_target_ligands)
        self.assertNotIn("J", sc.unassigned_target_ligands_reasons)

        # Check unassigned models
        # OXY: not identical to anything in the model
        mdl_oxy = sc.model.FindResidue("L_OXY", 1)
        self.assertEqual(sc.unassigned_model_ligands_reasons["L_OXY"][1], "identity")
        self.assertTrue("L_OXY" not in sc.score)
        # NA: not in contact with target
        mdl_na = sc.model.FindResidue("L_NA", 1)
        self.assertEqual(sc.unassigned_model_ligands_reasons["L_NA"][1], "no_contact")
        self.assertTrue("L_NA" not in sc.score)
 
        # MG in L_MG_2 has stupid coordinates and is not assigned
        mdl_mg_2 = sc.model.FindResidue("L_MG_2", 1)
        self.assertEqual(sc.unassigned_model_ligands_reasons["L_MG_2"][1], "stoichiometry")
        self.assertTrue("L_MG_2" not in sc.score)

        self.assertNotIn("L_MG_0", sc.unassigned_model_ligands_reasons)
        # CMO: disconnected
        mdl_cmo1 = sc.model.FindResidue("L_CMO", 1)
        self.assertEqual(sc.unassigned_model_ligands_reasons["L_CMO"][1], "disconnected")

        # Should work with rmsd_assignment too
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, trg_lig,
                                                full_bs_search=True)

        self.assertDictEqual(sc.unassigned_model_ligands_reasons, {
            'L_ZN': {1: 'model_binding_site'},
            'L_NA': {1: 'target_binding_site'},
            'L_OXY': {1: 'identity'},
            'L_MG_2': {1: 'stoichiometry'},
            "L_CMO": {1: 'disconnected'}
        })
        self.assertDictEqual(sc.unassigned_target_ligands_reasons, {
            'G': {1: 'identity'},
            'H': {1: 'model_binding_site'},
            'J': {1: 'stoichiometry'},
            'K': {1: 'identity'},
            'L_NA': {1: 'target_binding_site'},
            "L_CMO": {1: 'disconnected'}
        })
        self.assertTrue("L_OXY" not in sc.score)

        # With missing ligands
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, [], trg_lig)
        self.assertEqual(sc.unassigned_target_ligands_reasons["E"][1], 'no_ligand')

        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, mdl_lig, [])
        self.assertEqual(sc.unassigned_model_ligands_reasons["L_2"][1], 'no_ligand')

        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, [], trg_lig)
        self.assertEqual(sc.unassigned_target_ligands_reasons["E"][1], 'no_ligand')

        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, [])
        self.assertEqual(sc.unassigned_model_ligands_reasons["L_2"][1], 'no_ligand')

        # However not everything must be missing
        with self.assertRaises(ValueError):
            sc = LigandScorer(mdl, trg, [], [])

        # Test with partial bs search (full_bs_search=True)
        # Here we expect L_MG_2 to be unassigned because of stoichiometry
        # rather than model_binding_site, as it no longer matters so far from
        # the binding site.
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, trg_lig,
                                               full_bs_search=True)
        self.assertEqual(sc.unassigned_model_ligands_reasons, {
            'L_ZN': {1: 'model_binding_site'},
            'L_NA': {1: 'target_binding_site'},
            'L_OXY': {1: 'identity'},
            'L_MG_2': {1: 'stoichiometry'},
            "L_CMO": {1: 'disconnected'}
        })
        self.assertEqual(sc.unassigned_target_ligands_reasons, {
            'G': {1: 'identity'},
            'H': {1: 'model_binding_site'},
            'J': {1: 'stoichiometry'},
            'K': {1: 'identity'},
            'L_NA': {1: 'target_binding_site'},
            "L_CMO": {1: 'disconnected'}
        })

    def test_cleanup_polymer_ent(self):

        trg, trg_lig = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                              extract_nonpoly=True)
        mdl, mdl_lig = scoring_base.MMCIFPrep(_GetTestfilePath("P84080_model_02.cif.gz"),
                                              extract_nonpoly=True)

        # check hydrogen cleanup
        trg_with_h = trg.Copy()
        N = trg_with_h.GetAtomCount()
        ed = trg_with_h.EditXCS()
        ed.InsertAtom(trg_with_h.residues[0], "H", geom.Vec3(), "H")
        self.assertEqual(trg_with_h.GetAtomCount(), N + 1)
        sc = LigandScorer(mdl, trg_with_h, mdl_lig, trg_lig)
        self.assertEqual(sc.target.GetAtomCount(),N)

        trg_with_d = trg.Copy()
        N = trg_with_d.GetAtomCount()
        ed = trg_with_d.EditXCS()
        ed.InsertAtom(trg_with_d.residues[0], "H", geom.Vec3(), "D")
        self.assertEqual(trg_with_d.GetAtomCount(), N + 1)
        sc = LigandScorer(mdl, trg_with_d, mdl_lig, trg_lig)
        self.assertEqual(sc.target.GetAtomCount(),N)

        mdl_with_h = mdl.Copy()
        N = mdl_with_h.GetAtomCount()
        ed = mdl_with_h.EditXCS()
        ed.InsertAtom(mdl_with_h.residues[0], "H", geom.Vec3(), "H")
        self.assertEqual(mdl_with_h.GetAtomCount(), N + 1)
        sc = LigandScorer(mdl_with_h, trg, mdl_lig, trg_lig)
        self.assertEqual(sc.model.GetAtomCount(),N)

        mdl_with_d = mdl.Copy()
        N = mdl_with_d.GetAtomCount()
        ed = mdl_with_d.EditXCS()
        ed.InsertAtom(mdl_with_d.residues[0], "H", geom.Vec3(), "D")
        self.assertEqual(mdl_with_d.GetAtomCount(), N + 1)
        sc = LigandScorer(mdl_with_d, trg, mdl_lig, trg_lig)
        self.assertEqual(sc.model.GetAtomCount(),N)

        # residue with no entry in the component dictionary
        trg_with_funny_compound = trg.Copy()
        N = trg_with_funny_compound.GetResidueCount()
        ed = trg_with_funny_compound.EditXCS()
        ed.AppendResidue(trg_with_funny_compound.chains[0], "funny_compound")
        self.assertEqual(trg_with_funny_compound.GetResidueCount(), N+1)
        sc = LigandScorer(mdl, trg_with_funny_compound, mdl_lig, trg_lig)
        self.assertEqual(sc.target.GetResidueCount(),N)
        self.assertEqual(sc.target_cleanup_log["cleaned_residues"]["no_clib"], ["A.181."])

        mdl_with_funny_compound = trg.Copy()
        N = mdl_with_funny_compound.GetResidueCount()
        ed = mdl_with_funny_compound.EditXCS()
        ed.AppendResidue(mdl_with_funny_compound.chains[0], "funny_compound")
        self.assertEqual(mdl_with_funny_compound.GetResidueCount(), N+1)
        sc = LigandScorer(mdl_with_funny_compound, trg, mdl_lig, trg_lig)
        self.assertEqual(sc.model.GetResidueCount(),N)
        self.assertEqual(sc.model_cleanup_log["cleaned_residues"]["no_clib"], ["A.181."])

        # residue which is not peptide linking or nucleotide linking
        trg_with_pot = trg.Copy()
        N = trg_with_pot.GetResidueCount()
        ed = trg_with_pot.EditXCS()
        ed.AppendResidue(trg_with_pot.chains[0], "POT")
        self.assertTrue(conop.GetDefaultLib().FindCompound("POT") is not None)
        self.assertEqual(trg_with_pot.GetResidueCount(), N+1)
        sc = LigandScorer(mdl, trg_with_pot, mdl_lig, trg_lig)
        self.assertEqual(sc.target.GetResidueCount(),N)
        self.assertEqual(sc.target_cleanup_log["cleaned_residues"]["not_linking"], ["A.181."])

        mdl_with_pot = mdl.Copy()
        N = mdl_with_pot.GetResidueCount()
        ed = mdl_with_pot.EditXCS()
        ed.AppendResidue(mdl_with_pot.chains[0], "POT")
        self.assertTrue(conop.GetDefaultLib().FindCompound("POT") is not None)
        self.assertEqual(mdl_with_pot.GetResidueCount(), N+1)
        sc = LigandScorer(mdl_with_pot, trg, mdl_lig, trg_lig)
        self.assertEqual(sc.model.GetResidueCount(),N)
        self.assertEqual(sc.model_cleanup_log["cleaned_residues"]["not_linking"], ["A.181."])

        # unknown atom
        trg_with_unk = trg.Copy()
        N = trg_with_unk.GetAtomCount()
        N_res = trg_with_unk.GetResidueCount()
        ed = trg_with_unk.EditXCS()
        ed.InsertAtom(trg_with_unk.residues[0], "yolo", geom.Vec3(), "C")
        self.assertEqual(trg_with_unk.GetAtomCount(), N + 1)
        self.assertEqual(trg_with_unk.GetResidueCount(), N_res)
        sc = LigandScorer(mdl, trg_with_unk, mdl_lig, trg_lig)
        self.assertEqual(sc.target.GetAtomCount(), N)
        self.assertEqual(sc.target.GetResidueCount(), N_res)
        self.assertEqual(sc.target_cleanup_log["cleaned_atoms"]["unknown_atoms"], ["A.2..yolo"])        

        mdl_with_unk = mdl.Copy()
        N = mdl_with_unk.GetAtomCount()
        N_res = mdl_with_unk.GetResidueCount()
        ed = mdl_with_unk.EditXCS()
        ed.InsertAtom(mdl_with_unk.residues[0], "yolo", geom.Vec3(), "C")
        self.assertEqual(mdl_with_unk.GetAtomCount(), N + 1)
        self.assertEqual(mdl_with_unk.GetResidueCount(), N_res)
        sc = LigandScorer(mdl_with_unk, trg, mdl_lig, trg_lig)
        self.assertEqual(sc.model.GetAtomCount(), N)
        self.assertEqual(sc.model.GetResidueCount(), N_res)
        self.assertEqual(sc.model_cleanup_log["cleaned_atoms"]["unknown_atoms"], ["A.2..yolo"]) 

    def test_seqres(self):
        # the tested structures are absolutely irrelevant for ligands
        # we're simply checking if the seqres info gets propagated to
        # the chain mapper
        mdl = _LoadPDB("1eud_mdl_partial-dimer.pdb")
        trg = _LoadPDB("1eud_ref.pdb")
        mdl_lig = []
        trg_lig = []
        fake_lig = mol.CreateEntity()
        ed = fake_lig.EditXCS()
        ch=ed.InsertChain("A")
        r=ed.AppendResidue(ch, "LIG")
        a=ed.InsertAtom(r, "A", geom.Vec3())
        trg_lig.append(fake_lig)
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, trg_lig)
        self.assertTrue(sc._chain_mapper.seqres is None)
        self.assertTrue(sc._chain_mapper.trg_seqres_mapping is None)
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, mdl_lig, trg_lig)
        self.assertTrue(sc._chain_mapper.seqres is None)
        self.assertTrue(sc._chain_mapper.trg_seqres_mapping is None)


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
        sc = ligand_scoring_scrmsd.SCRMSDScorer(mdl, trg, mdl_lig, trg_lig,
                                                seqres=seqres,
                                                trg_seqres_mapping=trg_seqres_mapping,
                                                resnum_alignments=True)
        self.assertFalse(sc._chain_mapper.seqres is None)
        self.assertFalse(sc._chain_mapper.trg_seqres_mapping is None)
        sc = ligand_scoring_lddtpli.LDDTPLIScorer(mdl, trg, mdl_lig, trg_lig,
                                                  seqres=seqres,
                                                  trg_seqres_mapping=trg_seqres_mapping,
                                                  resnum_alignments=True)
        self.assertFalse(sc._chain_mapper.seqres is None)
        self.assertFalse(sc._chain_mapper.trg_seqres_mapping is None)


if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_ligand_scoring.py tests.')