import unittest, os, sys
from functools import lru_cache

import numpy as np

import ost
from ost import io, mol, geom
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.ligand_scoring import *
    from ost.mol.alg import ligand_scoring
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


class TestLigandScoring(unittest.TestCase):

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

        graph = ligand_scoring._ResidueToGraph(mdl_lig.residues[0])
        self.assertEqual(len(graph.edges), 34)
        self.assertEqual(len(graph.nodes), 32)
        # Check an arbitrary node
        self.assertEqual([a for a in graph.adj["14"].keys()], ["13", "29"])

        graph = ligand_scoring._ResidueToGraph(mdl_lig.residues[0], by_atom_index=True)
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

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1)
        self.assertEqual(len(sym), 72)

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1, by_atom_index=True)
        self.assertEqual(len(sym), 72)

        # Test that we can match ions read from SDF
        sdf_lig = _LoadEntity("1r8q_ligand_0.sdf")
        sym = ligand_scoring._ComputeSymmetries(trg_mg1, sdf_lig.residues[0], by_atom_index=True)
        self.assertEqual(len(sym), 1)

        # Test that it works with views and only consider atoms in the view
        # Skip PA, PB and O[1-3]A and O[1-3]B in target and model
        # We assume atom index are fixed and won't change
        trg_g3d1_sub_ent = trg_g3d1.Select("aindex>6019")
        trg_g3d1_sub = trg_g3d1_sub_ent.residues[0]
        mdl_g3d_sub_ent = mdl_g3d.Select("aindex>1447")
        mdl_g3d_sub = mdl_g3d_sub_ent.residues[0]

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub)
        self.assertEqual(len(sym), 6)

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub, by_atom_index=True)
        self.assertEqual(len(sym), 6)

        # Substructure matches
        sym = ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1_sub, substructure_match=True)
        self.assertEqual(len(sym), 6)

        # Missing atoms only allowed in target, not in model
        with self.assertRaises(NoSymmetryError):
            ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1, substructure_match=True)

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

        rmsd = SCRMSD(mdl_g3d, trg_g3d1)
        self.assertAlmostEqual(rmsd, 2.21341e-06, 10)
        rmsd = SCRMSD(mdl_g3d, trg_g3d2)
        self.assertAlmostEqual(rmsd, 61.21325, 4)

        # Ensure we raise a NoSymmetryError if the ligand is wrong
        with self.assertRaises(NoSymmetryError):
            SCRMSD(mdl_g3d, trg_mg1)
        with self.assertRaises(NoSymmetryError):
            SCRMSD(mdl_g3d, trg_afb1)

        # Assert that transform works
        trans = geom.Mat4(-0.999256, 0.00788487, -0.0377333, -15.4397,
                          0.0380652, 0.0473315, -0.998154, 29.9477,
                          -0.00608426, -0.998848, -0.0475963, 28.8251,
                          0, 0, 0, 1)
        rmsd = SCRMSD(mdl_g3d, trg_g3d2, transformation=trans)
        self.assertAlmostEqual(rmsd, 0.293972, 5)

        # Assert that substructure matches work
        trg_g3d1_sub = trg_g3d1.Select("aindex>6019").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        # mdl_g3d_sub = mdl_g3d.Select("aindex>1447").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        with self.assertRaises(NoSymmetryError):
            SCRMSD(mdl_g3d, trg_g3d1_sub)  # no full match

        # But partial match is OK
        rmsd = SCRMSD(mdl_g3d, trg_g3d1_sub, substructure_match=True)
        self.assertAlmostEqual(rmsd, 2.2376232209353475e-06, 8)

        # Ensure it doesn't work the other way around - ie incomplete model is invalid
        with self.assertRaises(NoSymmetryError):
            SCRMSD(trg_g3d1_sub, mdl_g3d)  # no full match

    def test__compute_scores(self):
        """Test that _compute_scores works.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")
        mdl_lig = io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))
        sc = LigandScorer(mdl, trg, [mdl_lig], None,
                          add_mdl_contacts = False,
                          lddt_pli_rmsd_binding_site = True)

        # Note: expect warning about Binding site of H.ZN1 not mapped to the model
        sc._compute_scores()

        # Check RMSD
        self.assertEqual(sc.rmsd_matrix.shape, (7, 1))
        np.testing.assert_almost_equal(sc.rmsd_matrix, np.array(
            [[np.nan],
            [0.04244993],
            [np.nan],
            [np.nan],
            [np.nan],
            [0.29399303],
            [np.nan]]), decimal=5)

        # Check lDDT-PLI
        self.assertEqual(sc.lddt_pli_matrix.shape, (7, 1))
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[0, 0]))
        self.assertAlmostEqual(sc.lddt_pli_matrix[1, 0], 0.99843, 5)
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[2, 0]))
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[3, 0]))
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[4, 0]))
        self.assertAlmostEqual(sc.lddt_pli_matrix[5, 0], 1.0)
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[6, 0]))

    def test_check_resnames(self):
        """Test that the check_resname argument works.

        When set to True, it should raise an error if any residue in the
        representation of the binding site in the model has a different
        name than in the reference. Here we manually modify a residue
        name to achieve that effect.
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
            sc = LigandScorer(mdl, trg, [mdl.FindResidue("I", 1)], [trg.FindResidue("I", 1)], check_resnames=True)
            sc._compute_scores()

        sc = LigandScorer(mdl, trg, [mdl.FindResidue("I", 1)], [trg.FindResidue("I", 1)], check_resnames=False)
        sc._compute_scores()

    def test__scores(self):
        """Test that the scores are computed correctly
        """
        # 4C0A has more ligands
        trg = _LoadMMCIF("1r8q.cif.gz")
        trg_4c0a = _LoadMMCIF("4c0a.cif.gz")
        sc = LigandScorer(trg, trg_4c0a, None, None, check_resnames=False,
                          add_mdl_contacts=False, lddt_pli_rmsd_binding_site = True)
        expected_keys = {"J", "F"}
        self.assertFalse(expected_keys.symmetric_difference(sc.rmsd.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.rmsd_details.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.lddt_pli.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.lddt_pli_details.keys()))
        # rmsd
        self.assertAlmostEqual(sc.rmsd["J"][mol.ResNum(1)], 0.8016608357429504, 5)
        self.assertAlmostEqual(sc.rmsd["F"][mol.ResNum(1)], 0.9286373257637024, 5)
        # rmsd_details
        self.assertEqual(sc.rmsd_details["J"][mol.ResNum(1)]["chain_mapping"], {'F': 'D', 'C': 'C'})
        self.assertEqual(len(sc.rmsd_details["J"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.rmsd_details["J"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.rmsd_details["J"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.rmsd_details["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.rmsd_details["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertEqual(sc.rmsd_details["F"][mol.ResNum(1)]["chain_mapping"], {'B': 'B', 'G': 'A'})
        self.assertEqual(len(sc.rmsd_details["F"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.rmsd_details["F"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.rmsd_details["F"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.rmsd_details["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.rmsd_details["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')

        # lddt_pli
        self.assertAlmostEqual(sc.lddt_pli["J"][mol.ResNum(1)], 0.9127105666156202, 5)
        self.assertAlmostEqual(sc.lddt_pli["F"][mol.ResNum(1)], 0.915929203539823, 5)
        # lddt_pli_details
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["lddt_pli_n_contacts"], 653)
        self.assertEqual(len(sc.lddt_pli_details["J"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["lddt_pli_n_contacts"], 678)
        self.assertEqual(len(sc.lddt_pli_details["F"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')

        # lddt_pli with added mdl contacts
        sc = LigandScorer(trg, trg_4c0a, None, None, check_resnames=False,
                          add_mdl_contacts=True)
        self.assertAlmostEqual(sc.lddt_pli["J"][mol.ResNum(1)], 0.8988340192043895, 5)
        self.assertAlmostEqual(sc.lddt_pli["F"][mol.ResNum(1)], 0.9039735099337749, 5)
        # lddt_pli_details
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["lddt_pli_n_contacts"], 729)
        self.assertEqual(len(sc.lddt_pli_details["J"][mol.ResNum(1)]["bs_ref_res"]), 63)
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["lddt_pli_n_contacts"], 755)
        self.assertEqual(len(sc.lddt_pli_details["F"][mol.ResNum(1)]["bs_ref_res"]), 62)
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')


    def test_rmsd_assignment(self):
        """Test that the RMSD-based assignment works.

        For RMSD, A: A results in a better chain mapping. However, C: A is a
        better global chain mapping from an lDDT perspective (and lDDT-PLI).
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        # By default, assignment differs between RMSD and lDDT-PLI in this
        # specific test case, so we can first ensure it does.
        # For now we skip as this is slow
        # sc = LigandScorer(mdl, trg, None, None)
        # assert sc.rmsd_details["L_2"][1]["target_ligand"] != sc.lddt_pli_details["L_2"][1]["target_ligand"]

        # RMSD assignment forces the same assignment
        sc = LigandScorer(mdl, trg, None, None, rmsd_assignment=True)
        self.assertEqual(sc.rmsd_details["L_2"][1]["target_ligand"], sc.lddt_pli_details["L_2"][1]["target_ligand"])

    def test_ignore_binding_site(self):
        """Test that we ignore non polymer stuff in the binding site.
         NOTE: we should consider changing this behavior in the future and take
         other ligands, peptides and short oligomers into account for superposition.
         When that's the case this test should be adapter
         """
        trg = _LoadMMCIF("1SSP.cif.gz")
        sc = LigandScorer(trg, trg, None, None)
        expected_bs_ref_res = ['C.GLY62', 'C.GLN63', 'C.ASP64', 'C.PRO65', 'C.TYR66', 'C.CYS76', 'C.PHE77', 'C.ASN123', 'C.HIS187']
        ost.PushVerbosityLevel(ost.LogLevel.Error)
        self.assertEqual([str(r) for r in sc.rmsd_details["D"][1]["bs_ref_res"]], expected_bs_ref_res)
        ost.PopVerbosityLevel()

    def test_unassigned_reasons(self):
        """Test reasons for being unassigned."""
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        def _AppendResidueWithBonds(ed, chain, old_res):
            new_res = ed.AppendResidue(chain, old_res.name)
            for old_atom in old_res.atoms:
                ed.InsertAtom(new_res, old_atom.name, old_atom.pos, old_atom.element,
                              old_atom.occupancy, old_atom.b_factor, old_atom.is_hetatom)
            for old_bond in old_atom.bonds:
                new_first = new_res.FindAtom(old_bond.first.name)
                new_second = new_res.FindAtom(old_bond.second.name)
                ed.Connect(new_first, new_second)
            return new_res

        # Add interesting ligands to model and target
        mdl_ed = mdl.EditXCS()
        trg_ed = trg.EditXCS()

        # Add ZN: representation in the model (chain missing in model)
        new_chain = mdl_ed.InsertChain("L_ZN")
        mdl_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = _AppendResidueWithBonds(mdl_ed, new_chain, trg.Select("rname=ZN").residues[0].handle)
        new_res.is_ligand = True

        # Add NA: not in contact with target
        new_chain = trg_ed.InsertChain("L_NA")
        trg_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = trg_ed.AppendResidue(new_chain, "NA")
        new_atom = trg_ed.InsertAtom(new_res, "NA", geom.Vec3(100, 100, 100), "NA")
        new_res.is_ligand = True
        new_chain = mdl_ed.InsertChain("L_NA")
        mdl_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = mdl_ed.AppendResidue(new_chain, "NA")
        new_atom = mdl_ed.InsertAtom(new_res, "NA", geom.Vec3(100, 100, 100), "NA")
        new_res.is_ligand = True

        # Add OXY: no symmetry/ not identical -
        new_chain = mdl_ed.InsertChain("L_OXY")
        mdl_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = mdl_ed.AppendResidue(new_chain, "OXY")
        new_atom1 = mdl_ed.InsertAtom(new_res, "O1", geom.Vec3(0, 0, 0), "O")
        new_atom2 = mdl_ed.InsertAtom(new_res, "O2", geom.Vec3(1, 1, 1), "O")
        mdl_ed.Connect(new_atom1, new_atom2)
        new_res.is_ligand = True

        # Add CMO: disconnected
        new_chain = mdl_ed.InsertChain("L_CMO")
        mdl_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = mdl_ed.AppendResidue(new_chain, "CMO")
        new_atom1 = mdl_ed.InsertAtom(new_res, "O", geom.Vec3(0, 0, 0), "O")
        new_atom2 = mdl_ed.InsertAtom(new_res, "C", geom.Vec3(1, 1, 1), "O")
        new_res.is_ligand = True
        new_chain = trg_ed.InsertChain("L_CMO")
        trg_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        new_res = trg_ed.AppendResidue(new_chain, "CMO")
        new_atom1 = trg_ed.InsertAtom(new_res, "O", geom.Vec3(0, 0, 0), "O")
        new_atom2 = trg_ed.InsertAtom(new_res, "C", geom.Vec3(1, 1, 1), "O")
        new_res.is_ligand = True

        # Add 3 MG in model: assignment/stoichiometry
        mg_pos = [
            geom.Vec3(3.871, 12.343, 44.485),
            geom.Vec3(3.871, 12.343, 44.485) + 1,
            geom.Vec3(3.871, 12.343, 44.485) + 100
        ]
        for i in range(3):
            new_chain = mdl_ed.InsertChain("L_MG_%d" % i)
            mdl_ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
            new_res = mdl_ed.AppendResidue(new_chain, "MG")
            new_atom = mdl_ed.InsertAtom(new_res, "MG", mg_pos[i], "MG")
            new_res.is_ligand = True

        mdl_ed.UpdateICS()
        trg_ed.UpdateICS()

        print("1")
        sc = LigandScorer(mdl, trg, None, None, unassigned=True,
                          full_bs_search=True,
                          add_mdl_contacts=False)
        print("done")
        # Check unassigned targets
        # NA: not in contact with target
        trg_na = sc.target.FindResidue("L_NA", 1)
        self.assertEqual(sc.unassigned_target_ligands["L_NA"][1], "binding_site")
        # ZN: no representation
        trg_zn = sc.target.FindResidue("H", 1)
        self.assertEqual(sc.unassigned_target_ligands["H"][1], "model_representation")
        # AFB: not identical to anything in the model
        trg_afb = sc.target.FindResidue("G", 1)
        self.assertEqual(sc.unassigned_target_ligands["G"][1], "identity")
        # F.G3D1: J.G3D1 assigned instead
        trg_fg3d1 = sc.target.FindResidue("F", 1)
        self.assertEqual(sc.unassigned_target_ligands["F"][1], "stoichiometry")
        # CMO: disconnected
        trg_cmo1 = sc.target.FindResidue("L_CMO", 1)
        self.assertEqual(sc.unassigned_target_ligands["L_CMO"][1], "disconnected")
        # J.G3D1: assigned to L_2.G3D1 => error
        trg_jg3d1 = sc.target.FindResidue("J", 1)
        with self.assertRaises(RuntimeError):
            sc._find_unassigned_target_ligand_reason(trg_jg3d1)
        self.assertNotIn("J", sc.unassigned_target_ligands)
        # Raises with an invalid ligand
        with self.assertRaises(ValueError):
            sc._find_unassigned_target_ligand_reason(sc.model_ligands[0])

        # Check unassigned models
        # OXY: not identical to anything in the model
        mdl_oxy = sc.model.FindResidue("L_OXY", 1)
        self.assertEqual(sc.unassigned_model_ligands["L_OXY"][1], "identity")
        self.assertIsNone(sc.lddt_pli["L_OXY"][1])
        # NA: not in contact with target
        mdl_na = sc.model.FindResidue("L_NA", 1)
        self.assertEqual(sc.unassigned_model_ligands["L_NA"][1], "binding_site")
        self.assertIsNone(sc.lddt_pli["L_NA"][1])
        # ZN: no representation
        mdl_zn = sc.model.FindResidue("L_ZN", 1)
        self.assertEqual(sc.unassigned_model_ligands["L_ZN"][1], "model_representation")
        self.assertIsNone(sc.lddt_pli["L_ZN"][1])
        # MG in L_MG_2 has stupid coordinates and is not assigned
        mdl_mg_2 = sc.model.FindResidue("L_MG_2", 1)
        self.assertEqual(sc.unassigned_model_ligands["L_MG_2"][1], "stoichiometry")
        self.assertIsNone(sc.lddt_pli["L_MG_2"][1])
        # MG in L_MG_0: assigned to I.MG1 => error
        mdl_mg_0 = sc.model.FindResidue("L_MG_0", 1)
        with self.assertRaises(RuntimeError):
            sc._find_unassigned_model_ligand_reason(mdl_mg_0)
        self.assertNotIn("L_MG_0", sc.unassigned_model_ligands)
        # CMO: disconnected
        mdl_cmo1 = sc.model.FindResidue("L_CMO", 1)
        self.assertEqual(sc.unassigned_model_ligands["L_CMO"][1], "disconnected")
        # Raises with an invalid ligand
        with self.assertRaises(ValueError):
            sc._find_unassigned_model_ligand_reason(sc.target_ligands[0])

        # Should work with rmsd_assignment too
        sc = LigandScorer(mdl, trg, None, None, unassigned=True,
                          rmsd_assignment=True, full_bs_search=True,
                          add_mdl_contacts=False,
                          lddt_pli_rmsd_binding_site=True)
        self.assertEqual(sc.unassigned_model_ligands, {
            'L_ZN': {1: 'model_representation'},
            'L_NA': {1: 'binding_site'},
            'L_OXY': {1: 'identity'},
            'L_MG_2': {1: 'stoichiometry'},
            "L_CMO": {1: 'disconnected'}
        })
        self.assertEqual(sc.unassigned_target_ligands, {
            'G': {1: 'identity'},
            'H': {1: 'model_representation'},
            'J': {1: 'stoichiometry'},
            'K': {1: 'identity'},
            'L_NA': {1: 'binding_site'},
            "L_CMO": {1: 'disconnected'}
        })
        self.assertIsNone(sc.lddt_pli["L_OXY"][1])

        # With missing ligands
        sc = LigandScorer(mdl.Select("cname=A"), trg, None, None)
        self.assertEqual(sc.unassigned_target_ligands["E"][1], 'no_ligand')

        sc = LigandScorer(mdl, trg.Select("cname=A"), None, None)
        self.assertEqual(sc.unassigned_model_ligands["L_2"][1], 'no_ligand')

        sc = LigandScorer(mdl.Select("cname=A"), trg, None, None,
                          unassigned=True, rmsd_assignment=True)
        self.assertEqual(sc.unassigned_target_ligands["E"][1], 'no_ligand')

        sc = LigandScorer(mdl, trg.Select("cname=A"), None, None,
                          unassigned=True, rmsd_assignment=True)
        self.assertEqual(sc.unassigned_model_ligands["L_2"][1], 'no_ligand')

        # However not everything must be missing
        with self.assertRaises(ValueError):
            sc = LigandScorer(mdl.Select("cname=A"), trg.Select("cname=A"), None, None,
                              unassigned=True, rmsd_assignment=True)

        # Test with partial bs search (full_bs_search=False)
        # Here we expect L_MG_2 to be unassigned because of model_representation
        # rather than stoichiometry, as it is so far from the binding site that
        # there is no longer a model binding site.
        sc = LigandScorer(mdl, trg, None, None, unassigned=True,
                          rmsd_assignment=True, full_bs_search=False)
        self.assertEqual(sc.unassigned_model_ligands, {
            'L_ZN': {1: 'model_representation'},
            'L_NA': {1: 'binding_site'},
            'L_OXY': {1: 'identity'},
            'L_MG_2': {1: 'model_representation'},
            "L_CMO": {1: 'disconnected'}
        })
        self.assertEqual(sc.unassigned_target_ligands, {
            'G': {1: 'identity'},
            'H': {1: 'model_representation'},
            'J': {1: 'stoichiometry'},
            'K': {1: 'identity'},
            'L_NA': {1: 'binding_site'},
            "L_CMO": {1: 'disconnected'}
        })


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

        # Substructure matches
        sc = LigandScorer(mdl.Select("protein=True"), trg.Select("protein=True"),
                          model_ligands=[mdl_g3d], target_ligands=[trg_g3d1_sub],
                          substructure_match=True)
        self.assertEqual(sc.rmsd_details["L_2"][1]["coverage"], 0.75)

    def test_6jyf(self):
        """6JYF initially caused issues in the CASP15-CAMEO/LIGATE paper where
         the ligand RET was wrongly assigned to short copies of OLA that float
          around and yielded higher scores.
          Here we test that this is resolved correctly."""
        mdl = _LoadPDB("6jyf_mdl.pdb")
        trg = _LoadMMCIF("6jyf_trg.cif")
        mdl_lig = _LoadEntity("6jyf_RET_pred.sdf")
        mdl_lig_full = _LoadEntity("6jyf_RET_pred_complete.sdf")

        # Problem is easily fixed by just prioritizing full coverage
        sc = LigandScorer(mdl, trg, model_ligands=[mdl_lig],
                          substructure_match=True)
        self.assertEqual(sc.rmsd_details['00001_'][1]["coverage"], 1.0)
        self.assertEqual(sc.rmsd_details['00001_'][1]["target_ligand"].name, "RET")
        self.assertAlmostEqual(sc.rmsd['00001_'][1], 15.56022, 4)
        self.assertTrue(np.array_equal(sc.coverage_matrix,
                              np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.3, 0.45, 0, 0, 0.55]]).transpose()))

        # We need to make sure that it also works if the match is partial.
        # For that we load the complete ligand incl. the O missing in target
        # with a coverage of around 95% only.
        sc = LigandScorer(mdl, trg, model_ligands=[mdl_lig_full],
                          substructure_match=True)
        self.assertTrue(sc.rmsd_details['00001_'][1]["coverage"] > 0.95)
        self.assertEqual(sc.rmsd_details['00001_'][1]["target_ligand"].name, "RET")
        self.assertAlmostEqual(sc.rmsd['00001_'][1], 15.56022, 4)

        # Next, we check that coverage_delta has an effect. With a large
        # delta of 0.5 we will assign to OLA which has a higher RMSD
        # but a coverage of 0.52 only.
        sc = LigandScorer(mdl, trg, model_ligands=[mdl_lig_full],
                          substructure_match=True,
                          coverage_delta=0.5)
        self.assertTrue(sc.rmsd_details['00001_'][1]["coverage"] > 0.5)
        self.assertEqual(sc.rmsd_details['00001_'][1]["target_ligand"].name, "OLA")
        self.assertAlmostEqual(sc.rmsd['00001_'][1], 6.13006878, 4)


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
        symmetries = ligand_scoring._ComputeSymmetries(mdl_lig[0], trg_lig[0], max_symmetries=100)
        assert len(symmetries) == 72
        with self.assertRaises(TooManySymmetriesError):
            ligand_scoring._ComputeSymmetries(mdl_lig[0], trg_lig[0], max_symmetries=10)

        # Check the unassignment
        sc = LigandScorer(mdl, trg, mdl_lig, trg_lig, max_symmetries=10)
        assert "L_2" not in sc.rmsd
        sc.unassigned_model_ligands["L_2"][1] == "symmetries"
        sc.unassigned_target_ligands["F"][1] == "symmetries"


    def test_no_binding_site(self):
        """
        Test the behavior when there's no binding site in proximity of
        the ligand. This test was introduced to identify some subtle issues
        with the ligand assignment that can cause it to enter an infinite
        loop when the data matrices are not filled properly.
        """
        trg = _LoadMMCIF("1r8q.cif.gz").Copy()
        mdl = trg.Copy()

        trg_zn = trg.FindResidue("H", 1)
        trg_g3d = trg.FindResidue("F", 1)


        # Move the zinc out of the reference binding site...
        ed = trg.EditXCS()
        ed.SetAtomPos(trg_zn.FindAtom("ZN"),
                      trg_zn.FindAtom("ZN").pos + geom.Vec3(6, 0, 0))
        # Remove some atoms from G3D to decrease coverage. This messed up
        # the assignment in the past.
        ed.DeleteAtom(trg_g3d.FindAtom("O6"))
        ed.UpdateICS()

        sc = LigandScorer(mdl, trg, target_ligands=[trg_zn, trg_g3d],
                          coverage_delta=0, substructure_match=True)

        self.assertTrue(np.isnan(sc.rmsd_matrix[0, 3]))
        self.assertEqual(sc.unassigned_target_ligands["H"][1], "binding_site")


    def test_no_lddt_pli_contact(self):
        """
        Test behaviour where a binding site has no lDDT-PLI contacts.

        We give two copies of the target ligand which have binding site atoms
        within radius=5A but no atoms at 4A. We set lddt_pli_radius=4 so that
        there are no contacts for the lDDT-PLI computation, and lDDT is None.

        We check that:
        - We don't get an lDDT-PLI assignment
        - Both target ligands are unassigned and have the
        - We get an RMSD assignment
        - The second copy of the target and model ligands ensure that the
          disambiguation code (which checks for the best lDDT-PLI when 2 RMSDs
          are identical) works in this case (where there is no lDDT-PLI to
          disambiguate the RMSDs).
        - We get lDDT-PLI = None with RMSD assignment
        """
        trg = _LoadPDB("T1118v1.pdb")
        trg_lig = _LoadEntity("T1118v1_001.sdf")
        mdl = _LoadPDB("T1118v1LG035_1.pdb")
        mdl_lig = _LoadEntity("T1118v1LG035_1_1_1.sdf")

        # Ensure it's unassigned in lDDT
        sc = LigandScorer(mdl, trg, [mdl_lig, mdl_lig], [trg_lig, trg_lig],
                          radius=5, lddt_pli_radius=4, rename_ligand_chain=True,
                          add_mdl_contacts=False)
        self.assertEqual(sc.lddt_pli,  {})
        self.assertEqual(sc.unassigned_target_ligands['00001_C'][1], "no_contact")
        self.assertEqual(sc.unassigned_target_ligands['00001_C_2'][1], "no_contact")
        self.assertEqual(sc.unassigned_model_ligands['00001_FE'][1], "no_contact")
        self.assertEqual(sc.unassigned_model_ligands['00001_FE_2'][1], "no_contact")
        # However RMSD should be assigned
        self.assertAlmostEqual(sc.rmsd['00001_FE'][1], 2.1703, 4)
        self.assertAlmostEqual(sc.rmsd['00001_FE_2'][1], 2.1703, 4)

        # Check with RMSD assignment
        sc2 = LigandScorer(mdl, trg, [mdl_lig, mdl_lig], [trg_lig, trg_lig],
                          radius=5, lddt_pli_radius=4, rename_ligand_chain=True,
                          rmsd_assignment=True)
        self.assertEqual(sc2.unassigned_target_ligands, {})
        self.assertEqual(sc2.unassigned_model_ligands, {})
        self.assertIsNone(sc2.lddt_pli['00001_FE'][1])
        self.assertIsNone(sc2.lddt_pli['00001_FE_2'][1])
        self.assertAlmostEqual(sc2.rmsd['00001_FE'][1], 2.1703, 4)
        self.assertAlmostEqual(sc2.rmsd['00001_FE_2'][1], 2.1703, 4)

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
        sc = LigandScorer(mdl, trg, add_mdl_contacts=False)
        self.assertAlmostEqual(sc.lddt_pli_matrix[0,0], 1.0, 5)

        # if added model contacts are considered, contributions from chain B are
        # perfectly reproduced but all contacts of ligand towards chain A are
        # added as penalty
        sc = LigandScorer(mdl, trg, add_mdl_contacts=True)

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
                
        self.assertAlmostEqual(sc.lddt_pli_matrix[0,0],
                               B_count/(A_count + B_count), 5)

        # Same as before but additionally we remove residue TRP.66 
        # from chain C in the target to test mapping magic...
        # Chain C is NOT in contact with the ligand but we only
        # add contacts from chain A as penalty that are mappable
        # to the closest chain with same sequence. That would be
        # chain C
        query = "cname=B,G or (cname=C and rnum!=66)"
        trg = mol.CreateEntityFromView(prot.Select(query), True)
        sc = LigandScorer(mdl, trg, add_mdl_contacts=True)

        TRP66_count = 0
        for a in lig.atoms:
            close_atoms = mdl.FindWithin(a.GetPos(), sc.lddt_pli_radius)
            for ca in close_atoms:
                cname = ca.GetChain().GetName()
                if cname == "A" and ca.GetResidue().GetNumber().GetNum() == 66:
                    TRP66_count += 1

        self.assertEqual(TRP66_count, 134)

        # remove TRP66_count from original penalty
        self.assertAlmostEqual(sc.lddt_pli_matrix[0,0],
                               B_count/(A_count + B_count - TRP66_count), 5)        

        # Move a random atom in the model from chain B towards the ligand center
        # chain B is also present in the target and interacts with the ligand,
        # but that atom would be far away and thus adds to the penalty. Since
        # the ligand is small enough, the number of added contacts should be
        # exactly the number of ligand atoms.
        mdl_ed = mdl.EditXCS()
        at = mdl.FindResidue("B", mol.ResNum(8)).FindAtom("NZ")
        mdl_ed.SetAtomPos(at, lig.geometric_center)
        sc = LigandScorer(mdl, trg, add_mdl_contacts=True)

        # compared to the last assertAlmostEqual, we add the number of ligand
        # atoms as additional penalties
        self.assertAlmostEqual(sc.lddt_pli_matrix[0,0],
                               B_count/(A_count + B_count - TRP66_count + \
                               lig.GetAtomCount()), 5)

if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_ligand_scoring.py tests.')
