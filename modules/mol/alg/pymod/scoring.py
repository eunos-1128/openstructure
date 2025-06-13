import os
from ost import mol
from ost import seq
from ost import io
from ost import conop
from ost import settings
from ost import geom
from ost import LogScript, LogInfo, LogDebug
from ost.mol.alg import lddt
from ost.mol.alg import qsscore
from ost.mol.alg import chain_mapping
from ost.mol.alg import stereochemistry
from ost.mol.alg import dockq
from ost.mol.alg.lddt import lDDTScorer
from ost.mol.alg.qsscore import QSScorer
from ost.mol.alg.contact_score import ContactScorer
from ost.mol.alg.contact_score import ContactEntity
from ost.mol.alg import GDT
from ost.mol.alg import Molck, MolckSettings
from ost import bindings
from ost.bindings import cadscore
from ost.bindings import tmtools
import numpy as np

def _GetAlignedResidues(aln, s1_ent, s2_ent):
    """ Yields aligned residues

    :param aln: The alignment with 2 sequences defining a residue-by-residue
                relationship.
    :type aln: :class:`ost.seq.AlignmentHandle`
    :param s1_ent: Structure representing first sequence in *aln*. 
                   One chain must be named after the first sequence and the
                   number of residues must match the number of non-gap
                   characters.
    :type s1_ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param s2_ent: Same for second sequence in *aln*.
    :type s2_ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    """
    s1_ch = s1_ent.FindChain(aln.GetSequence(0).GetName())
    s2_ch = s2_ent.FindChain(aln.GetSequence(1).GetName())

    if not s1_ch.IsValid():
        raise RuntimeError("s1_ent lacks required chain in _GetAlignedResidues")

    if not s2_ch.IsValid():
        raise RuntimeError("s2_ent lacks required chain in _GetAlignedResidues")

    if len(aln.GetSequence(0).GetGaplessString()) != s1_ch.GetResidueCount():
        raise RuntimeError("aln/chain mismatch in _GetAlignedResidues")
    if len(aln.GetSequence(1).GetGaplessString()) != s2_ch.GetResidueCount():
        raise RuntimeError("aln/chain mismatch in _GetAlignedResidues")

    s1_res = s1_ch.residues
    s2_res = s2_ch.residues

    s1_res_idx = 0
    s2_res_idx = 0

    for col in aln:
        if col[0] != '-' and col[1] != '-':
            yield (s1_res[s1_res_idx], s2_res[s2_res_idx])
        if col[0] != '-':
            s1_res_idx += 1
        if col[1] != '-':
            s2_res_idx += 1

class lDDTBSScorer:
    """Scorer specific for a reference/model pair

    Finds best possible binding site representation of reference in model given
    LDDT score. Uses :class:`ost.mol.alg.chain_mapping.ChainMapper` to deal with
    chain mapping.

    :param reference: Reference structure
    :type reference: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param model: Model structure
    :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param residue_number_alignment: Passed to ChainMapper constructor
    :type residue_number_alignment: :class:`bool`
    """
    def __init__(self, reference, model,
                 residue_number_alignment=False):
        self.chain_mapper = chain_mapping.ChainMapper(reference,
            resnum_alignments=residue_number_alignment)
        self.ref = self.chain_mapper.target
        self.mdl = model

    def ScoreBS(self, ligand, radius = 4.0, lddt_radius=10.0):
        """Computes binding site LDDT score given *ligand*. Best possible
        binding site representation is selected by LDDT but other scores such as
        CA based RMSD and GDT are computed too and returned.

        :param ligand: Defines the scored binding site, i.e. provides positions
                       to perform proximity search
        :type ligand: r'((Residue)|(Chain)|(Entity))((View)|(Handle))'
        :param radius: Reference residues with any atom position within *radius*
                       of *ligand* consitute the scored binding site
        :type radius: :class:`float`
        :param lddt_radius: Passed as *inclusion_radius* to
                            :class:`ost.mol.alg.lddt.lDDTScorer`
        :type lddt_radius: :class:`float`
        :returns: Object of type :class:`ost.mol.alg.chain_mapping.ReprResult`
                  containing all atom LDDT score and mapping information.
                  None if no representation could be found.
        """

        # create view of reference binding site
        ref_residues_hashes = set() # helper to keep track of added residues
        for ligand_at in ligand.atoms:
            close_atoms = self.ref.FindWithin(ligand_at.GetPos(), radius)
            for close_at in close_atoms:
                ref_res = close_at.GetResidue()
                h = ref_res.handle.GetHashCode()
                if h not in ref_residues_hashes:
                    ref_residues_hashes.add(h)

        # reason for doing that separately is to guarantee same ordering of
        # residues as in underlying entity. (Reorder by ResNum seems only
        # available on ChainHandles)
        ref_bs = self.ref.CreateEmptyView()
        for ch in self.ref.chains:
            for r in ch.residues:
                if r.handle.GetHashCode() in ref_residues_hashes:
                    ref_bs.AddResidue(r, mol.ViewAddFlag.INCLUDE_ALL)

        # gogogo
        bs_repr = self.chain_mapper.GetRepr(ref_bs, self.mdl,
                                            inclusion_radius = lddt_radius)
        if len(bs_repr) >= 1:
            return bs_repr[0]
        else:
            return None


class Scorer:
    """ Helper class to access the various scores available from ost.mol.alg

    Deals with structure cleanup, chain mapping, interface identification etc.
    Intermediate results are available as attributes.

    :param model: Model structure - a deep copy is available as :attr:`~model`.
                  Additionally, :func:`ost.mol.alg.Molck` using *molck_settings*
                  is applied.
    :type model: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param target: Target structure - a deep copy is available as :attr:`~target`.
                  Additionally, :func:`ost.mol.alg.Molck` using *molck_settings*
                  is applied.
    :type target: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param resnum_alignments: Whether alignments between chemically equivalent
                              chains in *model* and *target* can be computed
                              based on residue numbers. This can be assumed in
                              benchmarking setups such as CAMEO/CASP.
    :type resnum_alignments: :class:`bool`
    :param molck_settings: Settings used for Molck on *model* and *target*, if
                           set to None, a default object is constructed by
                           setting everything except rm_zero_occ_atoms and
                           colored to True in
                           :class:`ost.mol.alg.MolckSettings` constructor.
    :type molck_settings: :class:`ost.mol.alg.MolckSettings`
    :param cad_score_exec: Explicit path to voronota-cadscore executable from
                           voronota installation from 
                           https://github.com/kliment-olechnovic/voronota. If
                           not given, voronota-cadscore must be in PATH if any
                           of the CAD score related attributes is requested.
    :type cad_score_exec: :class:`str`
    :param custom_mapping: Provide custom chain mapping between *model* and
                           *target*. Dictionary with target chain names as key
                           and model chain names as value.
                           :attr:`~mapping` is constructed from this.
    :type custom_mapping: :class:`dict`
    :param custom_rigid_mapping: Provide custom chain mapping between *model*
                                 and *target*. Dictionary with target chain
                                 names as key and model chain names as value.
                                 :attr:`~rigid_mapping` is constructed from
                                 this.
    :type custom_rigid_mapping: :class:`dict`
    :param usalign_exec: Explicit path to USalign executable used to compute
                         TM-score. If not given, TM-score will be computed
                         with OpenStructure internal copy of USalign code.
    :type usalign_exec: :class:`str`
    :param lddt_no_stereochecks: Whether to compute LDDT without stereochemistry
                                checks
    :type lddt_no_stereochecks: :class:`bool`
    :param n_max_naive: Parameter for chain mapping. If the number of possible
                        mappings is <= *n_max_naive*, the full
                        mapping solution space is enumerated to find the
                        the optimum. A heuristic is used otherwise. The default
                        of 40320 corresponds to an octamer (8! = 40320).
                        A structure with stoichiometry A6B2 would be
                        6!*2! = 1440 etc.
    :type n_max_naive: :class:`int`
    :param oum: Override USalign Mapping. Inject rigid_mapping of
                :class:`Scorer` object into USalign to compute TM-score.
                Experimental feature with limitations.
    :type oum: :class:`bool`
    :param min_pep_length: Relevant parameter if short peptides are involved in
                           scoring. Minimum peptide length for a chain in the
                           target structure to be considered in chain mapping.
                           The chain mapping algorithm first performs an all vs.
                           all pairwise sequence alignment to identify \"equal\"
                           chains within the target structure. We go for simple
                           sequence identity there. Short sequences can be
                           problematic as they may produce high sequence identity
                           alignments by pure chance.
    :type min_pep_length: :class:`int`
    :param min_nuc_length: Relevant parameter if short nucleotides are involved
                           in scoring. Minimum nucleotide length for a chain in
                           the target structure to be considered in chain
                           mapping. The chain mapping algorithm first performs
                           an all vs. all pairwise sequence alignment to
                           identify \"equal\" chains within the target
                           structure. We go for simple sequence identity there.
                           Short sequences can be problematic as they may
                           produce high sequence identity alignments by pure
                           chance.
    :type min_nuc_length: :class:`int`
    :param lddt_add_mdl_contacts: LDDT specific flag. Only using contacts in
                                  LDDT that are within a certain distance
                                  threshold in the target does not penalize
                                  for added model contacts. If set to True, this
                                  flag will also consider target contacts
                                  that are within the specified distance
                                  threshold in the model but not necessarily in
                                  the target. No contact will be added if the
                                  respective atom pair is not resolved in the
                                  target.
    :type lddt_add_mdl_contacts: :class:`bool`
    :param dockq_capri_peptide: Flag that changes two things in the way DockQ
                                and its underlying scores are computed which is
                                proposed by the CAPRI community when scoring
                                peptides (PMID: 31886916).
                                ONE: Two residues are considered in contact if 
                                any of their atoms is within 5A. This is
                                relevant for fnat and fnonat scores. CAPRI
                                suggests to lower this threshold to 4A for
                                protein-peptide interactions.
                                TWO: irmsd is computed on interface residues. A
                                residue is defined as interface residue if any
                                of its atoms is within 10A of another chain.
                                CAPRI suggests to lower the default of 10A to
                                8A in combination with only considering CB atoms
                                for protein-peptide interactions.
                                This flag has no influence on patch_dockq
                                scores.
    :type dockq_capri_peptide: :class:`bool`
    :param lddt_symmetry_settings: Passed as *symmetry_settings* parameter to
                                   LDDT scorer. Default: None
    :type lddt_symmetry_settings: :class:`ost.mol.alg.lddt.SymmetrySettings`
    :param lddt_inclusion_radius: LDDT inclusion radius.
    :param pep_seqid_thr: Parameter that affects identification of identical
                          chains in target - see 
                          :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :type pep_seqid_thr: :class:`float`
    :param nuc_seqid_thr: Parameter that affects identification of identical
                          chains in target - see 
                          :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :type nuc_seqid_thr: :class:`float`
    :param mdl_map_pep_seqid_thr: Parameter that affects mapping of model chains
                                  to target chains - see 
                                  :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :type mdl_map_pep_seqid_thr: :class:`float`
    :param mdl_map_nuc_seqid_thr: Parameter that affects mapping of model chains
                                  to target chains - see 
                                  :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :type mdl_map_nuc_seqid_thr: :class:`float`
    :param seqres: Parameter that affects identification of identical chains in
                   target - see :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :type seqres: :class:`ost.seq.SequenceList`
    :param trg_seqres_mapping: Parameter that affects identification of identical
                               chains in target - see 
                               :class:`ost.mol.alg.chain_mapping.ChainMapper`
    :type trg_seqres_mapping: :class:`dict`
    """
    def __init__(self, model, target, resnum_alignments=False,
                 molck_settings = None, cad_score_exec = None,
                 custom_mapping=None, custom_rigid_mapping=None,
                 usalign_exec = None, lddt_no_stereochecks=False,
                 n_max_naive=40320, oum=False, min_pep_length = 6,
                 min_nuc_length = 4, lddt_add_mdl_contacts=False,
                 dockq_capri_peptide=False, lddt_symmetry_settings = None,
                 lddt_inclusion_radius = 15.0,
                 pep_seqid_thr = 95., nuc_seqid_thr = 95.,
                 mdl_map_pep_seqid_thr = 0.,
                 mdl_map_nuc_seqid_thr = 0.,
                 seqres = None,
                 trg_seqres_mapping = None):

        self._target_orig = target
        self._model_orig = model

        # lazily computed versions of target_orig and model_orig
        self._pepnuc_target = None
        self._pepnuc_model = None

        if isinstance(self._model_orig, mol.EntityView):
            self._model = mol.CreateEntityFromView(self._model_orig, False)
        else:
            self._model = self._model_orig.Copy()

        if isinstance(self._target_orig, mol.EntityView):
            self._target = mol.CreateEntityFromView(self._target_orig, False)
        else:
            self._target = self._target_orig.Copy()

        if molck_settings is None:
            molck_settings = MolckSettings(rm_unk_atoms=True,
                                           rm_non_std=False,
                                           rm_hyd_atoms=True,
                                           rm_oxt_atoms=True,
                                           rm_zero_occ_atoms=False,
                                           colored=False,
                                           map_nonstd_res=True,
                                           assign_elem=True)
        LogScript("Cleaning up input structures")
        Molck(self._model, conop.GetDefaultLib(), molck_settings)
        Molck(self._target, conop.GetDefaultLib(), molck_settings)

        if resnum_alignments:
            # If we're dealing with resnum alignments, we ensure that
            # consecutive peptide and nucleotide residues are connected based
            # on residue number information. The conop.Processor only connects
            # these things if the bonds are actually feasible which can lead to
            # weird behaviour in stereochemistry checks. Let's say N and C are
            # too close, it's reported as a clash. If they're too far apart,
            # they're not reported at all. If we're not dealing with resnum
            # alignments, we're out of luck as we have no direct residue
            # connectivity information from residue numbers
            self._resnum_connect(self._model)
            self._resnum_connect(self._target)

        self._model = self._model.Select("peptide=True or nucleotide=True")
        self._target = self._target.Select("peptide=True or nucleotide=True")

        # catch models which have empty chain names
        for ch in self._model.chains:
            if ch.GetName().strip() == "":
                raise RuntimeError("Model chains must have valid chain names")
            if ch.GetName().strip() == "'" or ch.GetName().strip() == '"':
                raise RuntimeError("Model chains cannot be named with quote signs (' or \"\")")
        
        # catch targets which have empty chain names
        for ch in self._target.chains:
            if ch.GetName().strip() == "":
                raise RuntimeError("Target chains must have valid chain names")
            if ch.GetName().strip() == "'" or ch.GetName().strip() == '"':
                raise RuntimeError("Target chains cannot be named with quote signs (' or \"\")")

        if resnum_alignments:
            # In case of resnum_alignments, we have some requirements on 
            # residue numbers in the chain mapping: 1) no ins codes 2) strictly
            # increasing residue numbers.
            for ch in self._model.chains:
                ins_codes = [r.GetNumber().GetInsCode() for r in ch.residues]
                if len(set(ins_codes)) != 1 or ins_codes[0] != '\0':
                    raise RuntimeError("Residue numbers in each model chain "
                                       "must not contain insertion codes if "
                                       "resnum_alignments are enabled")
                nums = [r.GetNumber().GetNum() for r in ch.residues]
                if not all(i < j for i, j in zip(nums, nums[1:])):
                    raise RuntimeError("Residue numbers in each model chain "
                                       "must be strictly increasing if "
                                       "resnum_alignments are enabled")

            for ch in self._target.chains:
                ins_codes = [r.GetNumber().GetInsCode() for r in ch.residues]
                if len(set(ins_codes)) != 1 or ins_codes[0] != '\0':
                    raise RuntimeError("Residue numbers in each target chain "
                                       "must not contain insertion codes if "
                                       "resnum_alignments are enabled")
                nums = [r.GetNumber().GetNum() for r in ch.residues]
                if not all(i < j for i, j in zip(nums, nums[1:])):
                    raise RuntimeError("Residue numbers in each target chain "
                                       "must be strictly increasing if "
                                       "resnum_alignments are enabled")

        if usalign_exec is not None:
            if not os.path.exists(usalign_exec):
                raise RuntimeError(f"USalign exec ({usalign_exec}) "
                                   f"not found")
            if not os.access(usalign_exec, os.X_OK):
                raise RuntimeError(f"USalign exec ({usalign_exec}) "
                                   f"is not executable")

        self.resnum_alignments = resnum_alignments
        self.cad_score_exec = cad_score_exec
        self.usalign_exec = usalign_exec
        self.lddt_no_stereochecks = lddt_no_stereochecks
        self.n_max_naive = n_max_naive
        self.oum = oum
        self.min_pep_length = min_pep_length
        self.min_nuc_length = min_nuc_length
        self.lddt_add_mdl_contacts = lddt_add_mdl_contacts
        self.dockq_capri_peptide = dockq_capri_peptide
        self.lddt_symmetry_settings = lddt_symmetry_settings
        self.lddt_inclusion_radius = lddt_inclusion_radius
        self.pep_seqid_thr = pep_seqid_thr
        self.nuc_seqid_thr = nuc_seqid_thr
        self.mdl_map_pep_seqid_thr = mdl_map_pep_seqid_thr
        self.mdl_map_nuc_seqid_thr = mdl_map_nuc_seqid_thr
        self.seqres = seqres
        self.trg_seqres_mapping = trg_seqres_mapping

        # lazily evaluated attributes
        self._stereochecked_model = None
        self._stereochecked_target = None
        self._model_clashes = None
        self._model_bad_bonds = None
        self._model_bad_angles = None
        self._target_clashes = None
        self._target_bad_bonds = None
        self._target_bad_angles = None
        self._trimmed_model = None
        self._chain_mapper = None
        self._mapping = None
        self._rigid_mapping = None
        self._model_interface_residues = None
        self._target_interface_residues = None
        self._aln = None
        self._stereochecked_aln = None
        self._pepnuc_aln = None
        self._trimmed_aln = None

        # lazily constructed scorer objects
        self._lddt_scorer = None
        self._bb_lddt_scorer = None
        self._qs_scorer = None
        self._contact_scorer = None
        self._trimmed_contact_scorer = None

        # lazily computed scores
        self._lddt = None
        self._local_lddt = None
        self._aa_local_lddt = None
        self._bb_lddt = None
        self._bb_local_lddt = None
        self._ilddt = None

        self._qs_global = None
        self._qs_best = None
        self._qs_target_interfaces = None
        self._qs_model_interfaces = None
        self._qs_interfaces = None
        self._per_interface_qs_global = None
        self._per_interface_qs_best = None

        self._contact_target_interfaces = None
        self._contact_model_interfaces = None
        self._native_contacts = None
        self._model_contacts = None
        self._trimmed_model_contacts = None
        self._ics_precision = None
        self._ics_recall = None
        self._ics = None
        self._per_interface_ics_precision = None
        self._per_interface_ics_recall = None
        self._per_interface_ics = None
        self._ips_precision = None
        self._ips_recall = None
        self._ips = None
        self._per_interface_ips_precision = None
        self._per_interface_ips_recall = None
        self._per_interface_ips = None

        # subset of contact scores that operate on trimmed model
        # i.e. no contacts from model residues that are not present in
        # target
        self._ics_trimmed = None
        self._ics_precision_trimmed = None
        self._ics_recall_trimmed = None
        self._per_interface_ics_precision_trimmed = None
        self._per_interface_ics_recall_trimmed = None
        self._per_interface_ics_trimmed = None
        self._ips_trimmed = None
        self._ips_precision_trimmed = None
        self._ips_recall_trimmed = None
        self._per_interface_ips_precision_trimmed = None
        self._per_interface_ips_recall_trimmed = None
        self._per_interface_ips_trimmed = None

        self._dockq_target_interfaces = None
        self._dockq_interfaces = None
        self._fnat = None
        self._fnonnat = None
        self._irmsd = None
        self._lrmsd = None
        self._nnat = None
        self._nmdl = None
        self._dockq_scores = None
        self._dockq_ave = None
        self._dockq_wave = None
        self._dockq_ave_full = None
        self._dockq_wave_full = None

        self._mapped_target_pos = None
        self._mapped_model_pos = None
        self._mapped_target_pos_full_bb = None
        self._mapped_model_pos_full_bb = None
        self._transformed_mapped_model_pos = None
        self._n_target_not_mapped = None
        self._transform = None

        self._rigid_mapped_target_pos = None
        self._rigid_mapped_model_pos = None
        self._rigid_mapped_target_pos_full_bb = None
        self._rigid_mapped_model_pos_full_bb = None
        self._rigid_transformed_mapped_model_pos = None
        self._rigid_n_target_not_mapped = None
        self._rigid_transform = None

        self._gdt_window_sizes = [7, 9, 12, 24, 48]
        self._gdt_05 = None
        self._gdt_1 = None
        self._gdt_2 = None
        self._gdt_4 = None
        self._gdt_8 = None
        self._gdtts = None
        self._gdtha = None
        self._rmsd = None

        self._cad_score = None
        self._local_cad_score = None

        self._patch_qs = None
        self._patch_dockq = None

        self._tm_score = None
        self._usalign_mapping = None

        if custom_mapping is not None:
            self._mapping = self._construct_custom_mapping(custom_mapping)

        if custom_rigid_mapping is not None:
            self._rigid_mapping = \
            self._construct_custom_mapping(custom_rigid_mapping)
        LogDebug("Scorer sucessfully initialized")

    @property
    def model(self):
        """ Model with Molck cleanup

        :type: :class:`ost.mol.EntityHandle`
        """
        return self._model

    @property
    def model_orig(self):
        """ The original model passed at object construction

        :type: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        """
        return self._model_orig

    @property
    def pepnuc_model(self):
        """ A selection of :attr:`~model_orig`

        Only contains peptide and nucleotide residues

        :type: :class:`ost.mol.EntityView`
        """
        if self._pepnuc_model is None:
            query = "peptide=true or nucleotide=true"
            self._pepnuc_model = self.model_orig.Select(query)
        return self._pepnuc_model

    @property
    def target(self):
        """ Target with Molck cleanup

        :type: :class:`ost.mol.EntityHandle`
        """
        return self._target

    @property
    def target_orig(self):
        """ The original target passed at object construction

        :type: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        """
        return self._target_orig

    @property
    def pepnuc_target(self):
        """ A selection of :attr:`~target_orig`

        Only contains peptide and nucleotide residues

        :type: :class:`ost.mol.EntityView`
        """
        if self._pepnuc_target is None:
            query = "peptide=true or nucleotide=true"
            self._pepnuc_target = self.target_orig.Select(query)
        return self._pepnuc_target

    @property
    def aln(self):
        """ Alignments of :attr:`~model`/:attr:`~target` chains

        Alignments for each pair of chains mapped in :attr:`~mapping`.
        First sequence is target sequence, second sequence the model sequence.

        :type: :class:`list` of :class:`ost.seq.AlignmentHandle`
        """
        if self._aln is None:
            self._compute_aln()
        return self._aln

    @property
    def stereochecked_aln(self):
        """ Stereochecked equivalent of :attr:`~aln`

        The alignments may differ, as stereochecks potentially remove residues

        :type: :class:`list` of :class:`ost.seq.AlignmentHandle`
        """
        if self._stereochecked_aln is None:
            self._compute_stereochecked_aln()
        return self._stereochecked_aln

    @property
    def pepnuc_aln(self):
        """ Alignments of :attr:`~model_orig`/:attr:`~target_orig` chains

        Selects for peptide and nucleotide residues before sequence
        extraction. Includes residues that would be removed by molck in
        structure preprocessing.

        :type: :class:`list` of :class:`ost.seq.AlignmentHandle`
        """
        if self._pepnuc_aln is None:
            self._compute_pepnuc_aln()
        return self._pepnuc_aln

    @property
    def trimmed_aln(self):
        """ Alignments of :attr:`~trimmed_model`/:attr:`~target` chains

        Alignments for each pair of chains mapped in :attr:`~mapping`.
        First sequence is target sequence, second sequence the model sequence.

        :type: :class:`list` of :class:`ost.seq.AlignmentHandle`
        """
        if self._trimmed_aln is None:
            self._trim_model()
        return self._trimmed_aln

    @property
    def stereochecked_model(self):
        """ View of :attr:`~model` that has stereochemistry checks applied

        First, a selection for peptide/nucleotide residues is performed,
        secondly peptide sidechains with stereochemical irregularities are
        removed (full residue if backbone atoms are involved). Irregularities
        are clashes or bond lengths/angles more than 12 standard deviations
        from expected values.

        :type: :class:`ost.mol.EntityView`
        """
        if self._stereochecked_model is None:
            self._do_stereochecks()
        return self._stereochecked_model

    @property
    def model_clashes(self):
        """ Clashing model atoms

        :type: :class:`list` of :class:`ost.mol.alg.stereochemistry.ClashInfo`
        """
        if self._model_clashes is None:
            self._do_stereochecks()
        return self._model_clashes

    @property
    def model_bad_bonds(self):
        """ Model bonds with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.BondViolationInfo`
        """
        if self._model_bad_bonds is None:
            self._do_stereochecks()
        return self._model_bad_bonds

    @property
    def model_bad_angles(self):
        """ Model angles with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.AngleViolationInfo`
        """
        if self._model_bad_angles is None:
            self._do_stereochecks()
        return self._model_bad_angles

    @property
    def stereochecked_target(self):
        """ Same as :attr:`~stereochecked_model` for :attr:`~target`

        :type: :class:`ost.mol.EntityView`
        """
        if self._stereochecked_target is None:
            self._do_stereochecks()
        return self._stereochecked_target

    @property
    def target_clashes(self):
        """ Clashing target atoms

        :type: :class:`list` of :class:`ost.mol.alg.stereochemistry.ClashInfo`
        """
        if self._target_clashes is None:
            self._do_stereochecks()
        return self._target_clashes

    @property
    def target_bad_bonds(self):
        """ Target bonds with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.BondViolationInfo`
        """
        if self._target_bad_bonds is None:
            self._do_stereochecks()
        return self._target_bad_bonds

    @property
    def target_bad_angles(self):
        """ Target angles with unexpected stereochemistry

        :type: :class:`list` of
               :class:`ost.mol.alg.stereochemistry.AngleViolationInfo`
        """
        if self._target_bad_angles is None:
            self._do_stereochecks()
        return self._target_bad_angles

    @property
    def trimmed_model(self):
        """ :attr:`~model` trimmed to target

        Removes residues that are not covered by :class:`target` given
        :attr:`~mapping`. In other words: no model residues without experimental
        evidence from :class:`target`. 

        :type: :class:`ost.mol.EntityView`
        """
        if self._trimmed_model is None:
            self._trim_model()
        return self._trimmed_model

    @property
    def chain_mapper(self):
        """ Chain mapper object for given :attr:`~target`

        :type: :class:`ost.mol.alg.chain_mapping.ChainMapper`
        """
        if self._chain_mapper is None:
            self._chain_mapper = chain_mapping.ChainMapper(self.target,
                                                           n_max_naive=1e9,
                                                           resnum_alignments=self.resnum_alignments,
                                                           min_pep_length=self.min_pep_length,
                                                           min_nuc_length=self.min_nuc_length,
                                                           pep_seqid_thr=self.pep_seqid_thr,
                                                           nuc_seqid_thr=self.nuc_seqid_thr,
                                                           mdl_map_pep_seqid_thr=self.mdl_map_pep_seqid_thr,
                                                           mdl_map_nuc_seqid_thr=self.mdl_map_nuc_seqid_thr,
                                                           seqres=self.seqres,
                                                           trg_seqres_mapping=self.trg_seqres_mapping)
        return self._chain_mapper

    @property
    def mapping(self):
        """ Full chain mapping result for :attr:`~target`/:attr:`~model`

        Computed with :func:`ost.mol.alg.ChainMapper.GetMapping`

        :type: :class:`ost.mol.alg.chain_mapping.MappingResult` 
        """
        if self._mapping is None:
            LogScript("Computing chain mapping")
            self._mapping = \
            self.chain_mapper.GetMapping(self.model,
                                         n_max_naive = self.n_max_naive)
        return self._mapping

    @property
    def rigid_mapping(self):
        """ Full chain mapping result for :attr:`~target`/:attr:`~model`

        Computed with :func:`ost.mol.alg.ChainMapper.GetRMSDMapping`

        :type: :class:`ost.mol.alg.chain_mapping.MappingResult` 
        """
        if self._rigid_mapping is None:
            LogScript("Computing rigid chain mapping")
            self._rigid_mapping = \
            self.chain_mapper.GetRMSDMapping(self.model)
        return self._rigid_mapping

    @property
    def model_interface_residues(self):
        """ Interface residues in :attr:`~model`

        Thats all residues having a contact with at least one residue from
        another chain (CB-CB distance <= 8A, CA in case of Glycine)

        :type: :class:`dict` with chain names as key and and :class:`list`
                with residue numbers of the respective interface residues.
        """
        if self._model_interface_residues is None:
            self._model_interface_residues = \
            self._get_interface_residues(self.model)
        return self._model_interface_residues

    @property
    def target_interface_residues(self):
        """ Same as :attr:`~model_interface_residues` for :attr:`~target`

        :type: :class:`dict` with chain names as key and and :class:`list`
                with residue numbers of the respective interface residues.
        """
        if self._target_interface_residues is None:
            self._target_interface_residues = \
            self._get_interface_residues(self.target)
        return self._target_interface_residues

    @property
    def lddt_scorer(self):
        """ LDDT scorer for :attr:`~target`/:attr:`~stereochecked_target`

        Depending on :attr:`~lddt_no_stereocheck` and
        :attr:`~lddt_symmetry_settings`.

        :type: :class:`ost.mol.alg.lddt.lDDTScorer`
        """
        if self._lddt_scorer is None:
            if self.lddt_no_stereochecks:
                self._lddt_scorer = lDDTScorer(self.target,
                                               symmetry_settings = self.lddt_symmetry_settings,
                                               inclusion_radius = self.lddt_inclusion_radius)
            else:
                self._lddt_scorer = lDDTScorer(self.stereochecked_target,
                                               symmetry_settings = self.lddt_symmetry_settings,
                                               inclusion_radius = self.lddt_inclusion_radius)
        return self._lddt_scorer

    @property
    def bb_lddt_scorer(self):
        """ LDDT scorer for :attr:`~target`, restricted to representative
        backbone atoms

        No stereochecks applied for bb only LDDT which considers CA atoms
        for peptides and C3' atoms for nucleotides.

        :type: :class:`ost.mol.alg.lddt.lDDTScorer`
        """
        if self._bb_lddt_scorer is None:
            self._bb_lddt_scorer = lDDTScorer(self.target, bb_only=True,
                                              symmetry_settings = self.lddt_symmetry_settings,
                                              inclusion_radius = self.lddt_inclusion_radius)
        return self._bb_lddt_scorer

    @property
    def qs_scorer(self):
        """ QS scorer constructed from :attr:`~mapping`

        The scorer object is constructed with default parameters and relates to
        :attr:`~model` and :attr:`~target` (no stereochecks).

        :type: :class:`ost.mol.alg.qsscore.QSScorer`
        """
        if self._qs_scorer is None:
            self._qs_scorer = QSScorer.FromMappingResult(self.mapping)
        return self._qs_scorer

    @property
    def contact_scorer(self):
        if self._contact_scorer is None:
            self._contact_scorer = ContactScorer.FromMappingResult(self.mapping)
        return self._contact_scorer
    
    @property
    def trimmed_contact_scorer(self):
        if self._trimmed_contact_scorer is None:
            self._trimmed_contact_scorer = ContactScorer(self.mapping.target,
                                                         self.mapping.chem_groups,
                                                         self.trimmed_model,
                                                         self.trimmed_aln)
        return self._trimmed_contact_scorer

    @property
    def lddt(self):
        """ Global LDDT score in range [0.0, 1.0]

        Computed based on :attr:`~stereochecked_model`. In case of oligomers,
        :attr:`~mapping` is used.

        :type: :class:`float`
        """
        if self._lddt is None:
            self._compute_lddt()
        return self._lddt
    
    @property
    def local_lddt(self):
        """ Per residue LDDT scores in range [0.0, 1.0]

        Computed based on :attr:`~stereochecked_model` but scores for all 
        residues in :attr:`~model` are reported. If a residue has been removed
        by stereochemistry checks, the respective score is set to 0.0. If a
        residue is not covered by the target or is in a chain skipped by the
        chain mapping procedure (happens for super short chains), the respective
        score is set to None. In case of oligomers, :attr:`~mapping` is used.

        :type: :class:`dict`
        """
        if self._local_lddt is None:
            self._compute_lddt()
        return self._local_lddt

    @property
    def aa_local_lddt(self):
        """ Per atom LDDT scores in range [0.0, 1.0]

        Computed based on :attr:`~stereochecked_model` but scores for all
        atoms in :attr:`~model` are reported. If an atom has been removed
        by stereochemistry checks, the respective score is set to 0.0. If an
        atom is not covered by the target or is in a chain skipped by the
        chain mapping procedure (happens for super short chains), the respective
        score is set to None. In case of oligomers, :attr:`~mapping` is used.

        :type: :class:`dict`
        """
        if self._aa_local_lddt is None:
            self._compute_lddt()
        return self._aa_local_lddt

    @property
    def bb_lddt(self):
        """ Global LDDT score restricted to representative backbone atoms in
        range [0.0, 1.0]

        Computed based on :attr:`~model` on representative backbone atoms only.
        This is CA for peptides and C3' for nucleotides. No stereochecks are
        performed. In case of oligomers, :attr:`~mapping` is used.

        :type: :class:`float`
        """
        if self._bb_lddt is None:
            self._compute_bb_lddt()
        return self._bb_lddt
    
    @property
    def bb_local_lddt(self):
        """ Per residue LDDT scores restricted to representative backbone atoms
        in range [0.0, 1.0]

        Computed based on :attr:`~model` on representative backbone atoms only.
        This is CA for peptides and C3' for nucleotides. No stereochecks are
        performed. If a residue is not covered by the target or is in a chain
        skipped by the chain mapping procedure (happens for super short
        chains), the respective score is set to None. In case of oligomers,
        :attr:`~mapping` is used.

        :type: :class:`dict`
        """
        if self._bb_local_lddt is None:
            self._compute_bb_lddt()
        return self._bb_local_lddt

    @property
    def ilddt(self):
        """ Global interface LDDT score in range [0.0, 1.0]

        This is LDDT only based on inter-chain contacts. Value is None if no
        such contacts are present. For example if we're dealing with a monomer.
        Computed based on :attr:`~stereochecked_model` and :attr:`~mapping` for
        chain mapping.

        :type: :class:`float`
        """
        if self._ilddt is None:
            # the whole None business kind of invalidates the idea of lazy
            # evaluation. The assumption is that this is called only once...
            self._compute_ilddt()
        return self._ilddt
    

    @property
    def qs_global(self):
        """  Global QS-score

        Computed based on :attr:`~model` using :attr:`~mapping`

        :type: :class:`float`
        """
        if self._qs_global is None:
            self._compute_qs()
        return self._qs_global

    @property
    def qs_best(self):
        """  Global QS-score - only computed on aligned residues

        Computed based on :attr:`~model` using :attr:`~mapping`. The QS-score
        computation only considers contacts between residues with a mapping
        between target and model. As a result, the score won't be lowered in
        case of additional chains/residues in any of the structures.

        :type: :class:`float`
        """
        if self._qs_best is None:
            self._compute_qs()
        return self._qs_best

    @property
    def qs_target_interfaces(self):
        """ Interfaces in :attr:`~target` with non-zero contribution to
        :attr:`~qs_global`/:attr:`~qs_best`

        Chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each:
               (trg_ch1, trg_ch2)
        """
        if self._qs_target_interfaces is None:
            self._qs_target_interfaces = self.qs_scorer.qsent1.interacting_chains
            self._qs_target_interfaces = \
            [(min(x[0],x[1]), max(x[0],x[1])) for x in self._qs_target_interfaces]
        return self._qs_target_interfaces

    @property
    def qs_model_interfaces(self):
        """ Interfaces in :attr:`~model` with non-zero contribution to
        :attr:`~qs_global`/:attr:`~qs_best`

        Chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each:
               (mdl_ch1, mdl_ch2)
        """
        if self._qs_model_interfaces is None:
            self._qs_model_interfaces = self.qs_scorer.qsent2.interacting_chains
            self._qs_model_interfaces = \
            [(min(x[0],x[1]), max(x[0],x[1])) for x in self._qs_model_interfaces]

        return self._qs_model_interfaces

    @property
    def qs_interfaces(self):
        """ Interfaces in :attr:`~qs_target_interfaces` that can be mapped
        to :attr:`~model`.

        Target chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 4 elements each:
               (trg_ch1, trg_ch2, mdl_ch1, mdl_ch2)
        """
        if self._qs_interfaces is None:
            self._qs_interfaces = list()
            flat_mapping = self.mapping.GetFlatMapping()
            for i in self.qs_target_interfaces:
                if i[0] in flat_mapping and i[1] in flat_mapping:
                    self._qs_interfaces.append((i[0], i[1],
                                                flat_mapping[i[0]],
                                                flat_mapping[i[1]]))
        return self._qs_interfaces
    
    @property
    def per_interface_qs_global(self):
        """ QS-score for each interface in :attr:`~qs_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_qs_global is None:
            self._compute_per_interface_qs_scores()
        return self._per_interface_qs_global
    
    @property
    def per_interface_qs_best(self):
        """ QS-score for each interface in :attr:`~qs_interfaces`

        Only computed on aligned residues

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_qs_best is None:
            self._compute_per_interface_qs_scores()
        return self._per_interface_qs_best
    
    @property
    def native_contacts(self):
        """ Native contacts

        A contact is a pair or residues from distinct chains that have
        a minimal heavy atom distance < 5A. Contacts are specified as
        :class:`tuple` with two strings in format:
        <cname>.<rnum>.<ins_code>

        :type: :class:`list` of :class:`tuple`
        """
        if self._native_contacts is None:
            self._native_contacts = self.contact_scorer.cent1.hr_contacts
        return self._native_contacts

    @property
    def model_contacts(self):
        """ Same for :attr:`~model`
        """
        if self._model_contacts is None:
            self._model_contacts = self.contact_scorer.cent2.hr_contacts
        return self._model_contacts

    @property
    def trimmed_model_contacts(self):
        """ Same for :attr:`~trimmed_model`
        """
        if self._trimmed_model_contacts is None:
            self._trimmed_model_contacts = self.trimmed_contact_scorer.cent2.hr_contacts
        return self._trimmed_model_contacts

    @property
    def contact_target_interfaces(self):
        """ Interfaces in :class:`target` which have at least one contact

        Contact as defined in :attr:`~native_contacts`,
        chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each
               (trg_ch1, trg_ch2)
        """
        if self._contact_target_interfaces is None:
            tmp = self.contact_scorer.cent1.interacting_chains
            tmp = [(min(x[0],x[1]), max(x[0],x[1])) for x in tmp]
            self._contact_target_interfaces = tmp
        return self._contact_target_interfaces

    @property
    def contact_model_interfaces(self):
        """ Interfaces in :class:`model` which have at least one contact

        Contact as defined in :attr:`~native_contacts`,
        chain names are lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each
               (mdl_ch1, mdl_ch2)
        """
        if self._contact_model_interfaces is None:
            tmp = self.contact_scorer.cent2.interacting_chains
            tmp = [(min(x[0],x[1]), max(x[0],x[1])) for x in tmp]
            self._contact_model_interfaces = tmp
        return self._contact_model_interfaces

    @property
    def ics_precision(self):
        """ Fraction of model contacts that are also present in target

        :type: :class:`float`
        """
        if self._ics_precision is None:
            self._compute_ics_scores()
        return self._ics_precision
    
    @property
    def ics_recall(self):
        """ Fraction of target contacts that are correctly reproduced in model

        :type: :class:`float`
        """
        if self._ics_recall is None:
            self._compute_ics_scores()
        return self._ics_recall

    @property
    def ics(self):
        """ ICS (Interface Contact Similarity) score

        Combination of :attr:`~ics_precision` and :attr:`~ics_recall`
        using the F1-measure

        :type: :class:`float`
        """
        if self._ics is None:
            self._compute_ics_scores()
        return self._ics

    @property
    def per_interface_ics_precision(self):
        """ Per-interface ICS precision

        :attr:`~ics_precision` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_precision is None:
            self._compute_ics_scores()
        return self._per_interface_ics_precision


    @property
    def per_interface_ics_recall(self):
        """ Per-interface ICS recall

        :attr:`~ics_recall` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_recall is None:
            self._compute_ics_scores()
        return self._per_interface_ics_recall

    @property
    def per_interface_ics(self):
        """ Per-interface ICS (Interface Contact Similarity) score

        :attr:`~ics` for each interface in 
        :attr:`~contact_target_interfaces`

        :type: :class:`float`
        """

        if self._per_interface_ics is None:
            self._compute_ics_scores()
        return self._per_interface_ics

    @property
    def ips_precision(self):
        """ Fraction of model interface residues that are also interface
        residues in target

        :type: :class:`float`
        """
        if self._ips_precision is None:
            self._compute_ips_scores()
        return self._ips_precision
    
    @property
    def ips_recall(self):
        """ Fraction of target interface residues that are also interface
        residues in model

        :type: :class:`float`
        """
        if self._ips_recall is None:
            self._compute_ips_scores()
        return self._ips_recall

    @property
    def ips(self):
        """ IPS (Interface Patch Similarity) score

        Jaccard coefficient of interface residues in target and their mapped
        counterparts in model

        :type: :class:`float`
        """
        if self._ips is None:
            self._compute_ips_scores()
        return self._ips

    @property
    def ics_trimmed(self):
        """ Same as :attr:`~ics` but with trimmed model

        Model is trimmed to residues which can me mapped to target in order
        to not penalize contacts in the model for which we have no experimental
        evidence.

        :type: :class:`float`
        """
        if self._ics_trimmed is None:
            self._compute_ics_scores_trimmed()
        return self._ics_trimmed

    @property
    def ics_precision_trimmed(self):
        """ Same as :attr:`~ics_precision` but with trimmed model

        Model is trimmed to residues which can me mapped to target in order
        to not penalize contacts in the model for which we have no experimental
        evidence.

        :type: :class:`float`
        """
        if self._ics_precision_trimmed is None:
            self._compute_ics_scores_trimmed()
        return self._ics_precision_trimmed

    @property
    def ics_recall_trimmed(self):
        """ Same as :attr:`~ics_recall` but with trimmed model

        Model is trimmed to residues which can me mapped to target in order
        to not penalize contacts in the model for which we have no experimental
        evidence.

        :type: :class:`float`
        """
        if self._ics_recall_trimmed is None:
            self._compute_ics_scores_trimmed()
        return self._ics_recall_trimmed

    @property
    def per_interface_ics_precision_trimmed(self):
        """ Same as :attr:`~per_interface_ics_precision` but with :attr:`~trimmed_model`

        :attr:`~ics_precision_trimmed` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_precision_trimmed is None:
            self._compute_ics_scores_trimmed()
        return self._per_interface_ics_precision_trimmed


    @property
    def per_interface_ics_recall_trimmed(self):
        """ Same as :attr:`~per_interface_ics_recall` but with :attr:`~trimmed_model`

        :attr:`~ics_recall_trimmed` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_recall_trimmed is None:
            self._compute_ics_scores_trimmed()
        return self._per_interface_ics_recall_trimmed

    @property
    def per_interface_ics_trimmed(self):
        """ Same as :attr:`~per_interface_ics` but with :attr:`~trimmed_model`

        :attr:`~ics` for each interface in 
        :attr:`~contact_target_interfaces`

        :type: :class:`float`
        """

        if self._per_interface_ics_trimmed is None:
            self._compute_ics_scores_trimmed()
        return self._per_interface_ics_trimmed

    @property
    def ips_trimmed(self):
        """ Same as :attr:`~ips` but with trimmed model

        Model is trimmed to residues which can me mapped to target in order
        to not penalize contacts in the model for which we have no experimental
        evidence.

        :type: :class:`float`
        """
        if self._ips_trimmed is None:
            self._compute_ips_scores_trimmed()
        return self._ips_trimmed

    @property
    def ips_precision_trimmed(self):
        """ Same as :attr:`~ips_precision` but with trimmed model

        Model is trimmed to residues which can me mapped to target in order
        to not penalize contacts in the model for which we have no experimental
        evidence.

        :type: :class:`float`
        """
        if self._ips_precision_trimmed is None:
            self._compute_ips_scores_trimmed()
        return self._ips_precision_trimmed

    @property
    def ips_recall_trimmed(self):
        """ Same as :attr:`~ips_recall` but with trimmed model

        Model is trimmed to residues which can me mapped to target in order
        to not penalize contacts in the model for which we have no experimental
        evidence.

        :type: :class:`float`
        """
        if self._ips_recall_trimmed is None:
            self._compute_ips_scores_trimmed()
        return self._ips_recall_trimmed

    @property
    def per_interface_ips_precision(self):
        """ Per-interface IPS precision

        :attr:`~ips_precision` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ips_precision is None:
            self._compute_ips_scores()
        return self._per_interface_ips_precision

    @property
    def per_interface_ips_recall(self):
        """ Per-interface IPS recall

        :attr:`~ips_recall` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ics_recall is None:
            self._compute_ips_scores()
        return self._per_interface_ips_recall

    @property
    def per_interface_ips(self):
        """ Per-interface IPS (Interface Patch Similarity) score

        :attr:`~ips` for each interface in 
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """

        if self._per_interface_ips is None:
            self._compute_ips_scores()
        return self._per_interface_ips

    @property
    def per_interface_ips_precision_trimmed(self):
        """ Same as :attr:`~per_interface_ips_precision` but with :attr:`~trimmed_model`

        :attr:`~ips_precision_trimmed` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ips_precision_trimmed is None:
            self._compute_ips_scores_trimmed()
        return self._per_interface_ips_precision_trimmed


    @property
    def per_interface_ips_recall_trimmed(self):
        """ Same as :attr:`~per_interface_ips_recall` but with :attr:`~trimmed_model`

        :attr:`~ics_recall_trimmed` for each interface in
        :attr:`~contact_target_interfaces`

        :type: :class:`list` of :class:`float`
        """
        if self._per_interface_ips_recall_trimmed is None:
            self._compute_ips_scores_trimmed()
        return self._per_interface_ips_recall_trimmed

    @property
    def per_interface_ips_trimmed(self):
        """ Same as :attr:`~per_interface_ips` but with :attr:`~trimmed_model`

        :attr:`~ics` for each interface in 
        :attr:`~contact_target_interfaces`

        :type: :class:`float`
        """

        if self._per_interface_ips_trimmed is None:
            self._compute_ips_scores_trimmed()
        return self._per_interface_ips_trimmed

    @property
    def dockq_target_interfaces(self):
        """ Interfaces in :attr:`~target` that are relevant for DockQ

        All interfaces in :attr:`~target` with non-zero contacts that are
        relevant for DockQ. Includes protein-protein, protein-nucleotide and
        nucleotide-nucleotide interfaces. Chain names for each interface are
        lexicographically sorted.

        :type: :class:`list` of :class:`tuple` with 2 elements each:
               (trg_ch1, trg_ch2)
        """
        if self._dockq_target_interfaces is None:
            
            # interacting chains are identified with ContactEntity
            contact_d = 5.0
            if self.dockq_capri_peptide:
                contact_d = 4.0
            cent = ContactEntity(self.target, contact_mode = "aa",
                                 contact_d = contact_d)

            # fetch lexicographically sorted interfaces
            interfaces = cent.interacting_chains
            interfaces = [(min(x[0],x[1]), max(x[0],x[1])) for x in interfaces]

            pep_seqs = set([s.GetName() for s in self.chain_mapper.polypep_seqs])
            nuc_seqs = set([s.GetName() for s in self.chain_mapper.polynuc_seqs])

            seqs = pep_seqs.union(nuc_seqs)
            self._dockq_target_interfaces = list()
            for interface in interfaces:
                if interface[0] in seqs and interface[1] in seqs:
                    self._dockq_target_interfaces.append(interface)

        return self._dockq_target_interfaces

    @property
    def dockq_interfaces(self):
        """ Interfaces in :attr:`~dockq_target_interfaces` that can be mapped
        to model

        Target chain names are lexicographically sorted

        :type: :class:`list` of :class:`tuple` with 4 elements each:
               (trg_ch1, trg_ch2, mdl_ch1, mdl_ch2)
        """
        if self._dockq_interfaces is None:
            self._dockq_interfaces = list()
            flat_mapping = self.mapping.GetFlatMapping()
            for i in self.dockq_target_interfaces:
                if i[0] in flat_mapping and i[1] in flat_mapping:
                    self._dockq_interfaces.append((i[0], i[1],
                                                   flat_mapping[i[0]],
                                                   flat_mapping[i[1]]))
        return self._dockq_interfaces
    
    @property
    def dockq_scores(self):
        """ DockQ scores for interfaces in :attr:`~dockq_interfaces` 

        :class:`list` of :class:`float`
        """
        if self._dockq_scores is None:
            self._compute_dockq_scores()
        return self._dockq_scores

    @property
    def fnat(self):
        """ fnat scores for interfaces in :attr:`~dockq_interfaces` 

        fnat: Fraction of native contacts that are also present in model

        :class:`list` of :class:`float`
        """
        if self._fnat is None:
            self._compute_dockq_scores()
        return self._fnat

    @property
    def nnat(self):
        """ N native contacts for interfaces in :attr:`~dockq_interfaces` 

        :class:`list` of :class:`int`
        """
        if self._nnat is None:
            self._compute_dockq_scores()
        return self._nnat

    @property
    def nmdl(self):
        """ N model contacts for interfaces in :attr:`~dockq_interfaces` 

        :class:`list` of :class:`int`
        """
        if self._nmdl is None:
            self._compute_dockq_scores()
        return self._nmdl

    @property
    def fnonnat(self):
        """ fnonnat scores for interfaces in :attr:`~dockq_interfaces` 

        fnat: Fraction of model contacts that are not present in target

        :class:`list` of :class:`float`
        """
        if self._fnonnat is None:
            self._compute_dockq_scores()
        return self._fnonnat

    @property
    def irmsd(self):
        """ irmsd scores for interfaces in :attr:`~dockq_interfaces` 

        irmsd: RMSD of interface (RMSD computed on backbone atoms) which
        consists of each residue that has at least one heavy atom within 10A of
        other chain. Backbone atoms for proteins: "CA","C","N","O", for
        nucleotides: "P", "OP1", "OP2", "O2'", "O3'", "O4'", "O5'", "C1'",
        "C2'", "C3'", "C4'", "C5'".

        :class:`list` of :class:`float`
        """
        if self._irmsd is None:
            self._compute_dockq_scores()
        return self._irmsd

    @property
    def lrmsd(self):
        """ lrmsd scores for interfaces in :attr:`~dockq_interfaces` 

        lrmsd: The two chains involved in the interface are superposed based on
        the receptor (rigid min RMSD superposition) and the ligand RMSD is
        reported. Receptor is the chain with more residues. Superposition and
        RMSD is computed on same backbone atoms as :attr:`~irmsd`.

        :class:`list` of :class:`float`
        """
        if self._lrmsd is None:
            self._compute_dockq_scores()
        return self._lrmsd
        
    @property
    def dockq_ave(self):
        """ Average of DockQ scores in :attr:`~dockq_scores`

        In its original implementation, DockQ only operates on single
        interfaces. Thus the requirement to combine scores for higher order
        oligomers.

        :type: :class:`float`
        """
        if self._dockq_ave is None:
            self._compute_dockq_scores()
        return self._dockq_ave
    
    @property
    def dockq_wave(self):
        """ Same as :attr:`~dockq_ave`, weighted by native contacts

        :type: :class:`float`
        """
        if self._dockq_wave is None:
            self._compute_dockq_scores()
        return self._dockq_wave
        
    @property
    def dockq_ave_full(self):
        """ Same as :attr:`~dockq_ave` but penalizing for missing interfaces

        Interfaces that are not covered in model are added as 0.0
        in average computation.

        :type: :class:`float`
        """
        if self._dockq_ave_full is None:
            self._compute_dockq_scores()
        return self._dockq_ave_full
    
    @property
    def dockq_wave_full(self):
        """ Same as :attr:`~dockq_ave_full`, but weighted

        Interfaces that are not covered in model are added as 0.0 in
        average computations and the respective weights are derived from
        number of contacts in respective target interface. 
        """
        if self._dockq_wave_full is None:
            self._compute_dockq_scores()
        return self._dockq_wave_full

    @property
    def mapped_target_pos(self):
        """ Mapped representative positions in target

        Thats CA positions for peptide residues and C3' positions for
        nucleotides. Has same length as :attr:`~mapped_model_pos` and mapping
        is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._mapped_target_pos is None:
            self._extract_mapped_pos()
        return self._mapped_target_pos

    @property
    def mapped_target_pos_full_bb(self):
        """ Mapped representative positions in target

        Thats the equivalent of :attr:`~mapped_target_pos` but containing more
        backbone atoms (N, CA, C for peptide residues and O5', C5', C4', C3', O3
        for nucleotide residues). mapping is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._mapped_target_pos_full_bb is None:
            self._extract_mapped_pos_full_bb()
        return self._mapped_target_pos_full_bb

    @property
    def mapped_model_pos(self):
        """ Mapped representative positions in model

        Thats CA positions for peptide residues and C3' positions for
        nucleotides. Has same length as :attr:`~mapped_target_pos` and mapping
        is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._mapped_model_pos is None:
            self._extract_mapped_pos()
        return self._mapped_model_pos

    @property
    def mapped_model_pos_full_bb(self):
        """ Mapped representative positions in model

        Thats the equivalent of :attr:`~mapped_model_pos` but containing more
        backbone atoms (N, CA, C for peptide residues and O5', C5', C4', C3', O3
        for nucleotide residues). mapping is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._mapped_model_pos_full_bb is None:
            self._extract_mapped_pos_full_bb()
        return self._mapped_model_pos_full_bb

    @property
    def transformed_mapped_model_pos(self):
        """ :attr:`~mapped_model_pos` with :attr:`~transform` applied

        :type: :class:`ost.geom.Vec3List`
        """
        if self._transformed_mapped_model_pos is None:
            self._transformed_mapped_model_pos = \
            geom.Vec3List(self.mapped_model_pos)
            self._transformed_mapped_model_pos.ApplyTransform(self.transform)
        return self._transformed_mapped_model_pos

    @property
    def n_target_not_mapped(self):
        """ Number of target residues which have no mapping to model

        :type: :class:`int`
        """
        if self._n_target_not_mapped is None:
            self._extract_mapped_pos()
        return self._n_target_not_mapped

    @property
    def transform(self):
        """ Transform: :attr:`~mapped_model_pos` onto :attr:`~mapped_target_pos`

        Computed using Kabsch minimal rmsd algorithm. If number of positions
        is too small (< 3), :attr:`~mapped_model_pos_full_bb` and
        :attr:`~mapped_target_pos_full_bb` are used.

        :type: :class:`ost.geom.Mat4`
        """
        if self._transform is None:
            if len(self.mapped_model_pos) < 3:
                if len(self.mapped_model_pos_full_bb) >=3:
                    res = mol.alg.SuperposeSVD(self.mapped_model_pos_full_bb,
                                               self.mapped_target_pos_full_bb)
                    self._transform = res.transformation
                else:
                    # there is really nothing we can do => set identity matrix
                    self._transform = geom.Mat4()
            else:
                res = mol.alg.SuperposeSVD(self.mapped_model_pos,
                                           self.mapped_target_pos)
                self._transform = res.transformation
        return self._transform

    @property
    def rigid_mapped_target_pos(self):
        """ Mapped representative positions in target

        Thats CA positions for peptide residues and C3' positions for
        nucleotides. Has same length as :attr:`~rigid_mapped_model_pos` and mapping
        is based on :attr:`~rigid_mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._rigid_mapped_target_pos is None:
            self._extract_rigid_mapped_pos()
        return self._rigid_mapped_target_pos

    @property
    def rigid_mapped_target_pos_full_bb(self):
        """ Mapped representative positions in target

        Thats the equivalent of :attr:`~rigid_mapped_target_pos` but containing
        more backbone atoms (N, CA, C for peptide residues and O5', C5', C4',
        C3', O3 for nucleotide residues). mapping is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._rigid_mapped_target_pos_full_bb is None:
            self._extract_rigid_mapped_pos_full_bb()
        return self._rigid_mapped_target_pos_full_bb

    @property
    def rigid_mapped_model_pos(self):
        """ Mapped representative positions in model

        Thats CA positions for peptide residues and C3' positions for
        nucleotides. Has same length as :attr:`~mapped_target_pos` and mapping
        is based on :attr:`~rigid_mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._rigid_mapped_model_pos is None:
            self._extract_rigid_mapped_pos()
        return self._rigid_mapped_model_pos

    @property
    def rigid_mapped_model_pos_full_bb(self):
        """ Mapped representative positions in model

        Thats the equivalent of :attr:`~rigid_mapped_model_pos` but containing
        more backbone atoms (N, CA, C for peptide residues and O5', C5', C4',
        C3', O3 for nucleotide residues). mapping is based on :attr:`~mapping`.

        :type: :class:`ost.geom.Vec3List`
        """
        if self._rigid_mapped_model_pos_full_bb is None:
            self._extract_rigid_mapped_pos_full_bb()
        return self._rigid_mapped_model_pos_full_bb

    @property
    def rigid_transformed_mapped_model_pos(self):
        """ :attr:`~rigid_mapped_model_pos` with :attr:`~rigid_transform` applied

        :type: :class:`ost.geom.Vec3List`
        """
        if self._rigid_transformed_mapped_model_pos is None:
            self._rigid_transformed_mapped_model_pos = \
            geom.Vec3List(self.rigid_mapped_model_pos)
            self._rigid_transformed_mapped_model_pos.ApplyTransform(self.rigid_transform)
        return self._rigid_transformed_mapped_model_pos

    @property
    def rigid_n_target_not_mapped(self):
        """ Number of target residues which have no rigid mapping to model

        :type: :class:`int`
        """
        if self._rigid_n_target_not_mapped is None:
            self._extract_rigid_mapped_pos()
        return self._rigid_n_target_not_mapped

    @property
    def rigid_transform(self):
        """ Transform: :attr:`~rigid_mapped_model_pos` onto :attr:`~rigid_mapped_target_pos`

        Computed using Kabsch minimal rmsd algorithm. If number of positions
        is too small (< 3), :attr:`~rigid_mapped_model_pos_full_bb` and
        :attr:`~rigid_mapped_target_pos_full_bb` are used.

        :type: :class:`ost.geom.Mat4`
        """
        if self._rigid_transform is None:
            if len(self.rigid_mapped_model_pos) < 3:
                if len(self.rigid_mapped_model_pos_full_bb) >= 3:
                    res = mol.alg.SuperposeSVD(self.rigid_mapped_model_pos_full_bb,
                                               self.rigid_mapped_target_pos_full_bb)
                    self._rigid_transform = res.transformation
                else:
                    # there is really nothing we can do => set identity matrix
                    self._rigid_transform = geom.Mat4()
            else:
                res = mol.alg.SuperposeSVD(self.rigid_mapped_model_pos,
                                           self.rigid_mapped_target_pos)
                self._rigid_transform = res.transformation
        return self._rigid_transform

    @property
    def gdt_05(self):
        """ Fraction CA (C3' for nucleotides) that can be superposed within 0.5A

        Uses :attr:`~rigid_mapped_model_pos` and :attr:`~rigid_mapped_target_pos`.
        Similar iterative algorithm as LGA tool

        :type: :class:`float` 
        """
        if self._gdt_05 is None:
            N = list()
            wsizes = self._gdt_window_sizes + [len(self.rigid_mapped_model_pos)]
            for window_size in wsizes:
                n = GDT(self.rigid_mapped_model_pos,
                        self.rigid_mapped_target_pos,
                        window_size, 1000, 0.5)[0]
                N.append(n)
            n = max(N)
            n_full = len(self.rigid_mapped_target_pos) + self.rigid_n_target_not_mapped
            if n_full > 0:
                self._gdt_05 = float(n) / n_full
            else:
                self._gdt_05 = 0.0
        return self._gdt_05

    @property
    def gdt_1(self):
        """ Fraction CA (C3' for nucleotides) that can be superposed within 1.0A

        Uses :attr:`~rigid_mapped_model_pos` and :attr:`~rigid_mapped_target_pos`.
        Similar iterative algorithm as LGA tool

        :type: :class:`float` 
        """
        if self._gdt_1 is None:
            N = list()
            wsizes = self._gdt_window_sizes + [len(self.rigid_mapped_model_pos)]
            for window_size in wsizes:
                n = GDT(self.rigid_mapped_model_pos,
                        self.rigid_mapped_target_pos,
                        window_size, 1000, 1.0)[0]
                N.append(n)
            n = max(N)
            n_full = len(self.rigid_mapped_target_pos) + self.rigid_n_target_not_mapped
            if n_full > 0:
                self._gdt_1 = float(n) / n_full
            else:
                self._gdt_1 = 0.0
        return self._gdt_1

    @property
    def gdt_2(self):
        """ Fraction CA (C3' for nucleotides) that can be superposed within 2.0A

        Uses :attr:`~rigid_mapped_model_pos` and :attr:`~rigid_mapped_target_pos`.
        Similar iterative algorithm as LGA tool


        :type: :class:`float` 
        """
        if self._gdt_2 is None:
            N = list()
            wsizes = self._gdt_window_sizes + [len(self.rigid_mapped_model_pos)]
            for window_size in wsizes:
                n = GDT(self.rigid_mapped_model_pos,
                        self.rigid_mapped_target_pos,
                        window_size, 1000, 2.0)[0]
                N.append(n)
            n = max(N)
            n_full = len(self.rigid_mapped_target_pos) + self.rigid_n_target_not_mapped
            if n_full > 0:
                self._gdt_2 = float(n) / n_full
            else:
                self._gdt_2 = 0.0
        return self._gdt_2

    @property
    def gdt_4(self):
        """ Fraction CA (C3' for nucleotides) that can be superposed within 4.0A

        Uses :attr:`~rigid_mapped_model_pos` and :attr:`~rigid_mapped_target_pos`.
        Similar iterative algorithm as LGA tool

        :type: :class:`float` 
        """
        if self._gdt_4 is None:
            N = list()
            wsizes = self._gdt_window_sizes + [len(self.rigid_mapped_model_pos)]
            for window_size in wsizes:
                n = GDT(self.rigid_mapped_model_pos,
                        self.rigid_mapped_target_pos,
                        window_size, 1000, 4.0)[0]
                N.append(n)
            n = max(N)
            n_full = len(self.rigid_mapped_target_pos) + self.rigid_n_target_not_mapped
            if n_full > 0:
                self._gdt_4 = float(n) / n_full
            else:
                self._gdt_4 = 0.0
        return self._gdt_4

    @property
    def gdt_8(self):
        """ Fraction CA (C3' for nucleotides) that can be superposed within 8.0A

        Similar iterative algorithm as LGA tool

        :type: :class:`float` 
        """
        if self._gdt_8 is None:
            N = list()
            wsizes = self._gdt_window_sizes + [len(self.rigid_mapped_model_pos)]
            for window_size in wsizes:
                n = GDT(self.rigid_mapped_model_pos,
                        self.rigid_mapped_target_pos,
                        window_size, 1000, 8.0)[0]
                N.append(n)
            n = max(N)
            n_full = len(self.rigid_mapped_target_pos) + self.rigid_n_target_not_mapped
            if n_full > 0:
                self._gdt_8 = float(n) / n_full
            else:
                self._gdt_8 = 0.0
        return self._gdt_8
    

    @property
    def gdtts(self):
        """ avg GDT with thresholds: 8.0A, 4.0A, 2.0A and 1.0A

        :type: :class:`float`
        """
        if self._gdtts is None:
            LogScript("Computing GDT-TS score")
            self._gdtts = (self.gdt_1 + self.gdt_2 + self.gdt_4 + self.gdt_8) / 4
        return self._gdtts

    @property
    def gdtha(self):
        """ avg GDT with thresholds: 4.0A, 2.0A, 1.0A and 0.5A

        :type: :class:`float`
        """
        if self._gdtha is None:
            LogScript("Computing GDT-HA score")
            self._gdtha = (self.gdt_05 + self.gdt_1 + self.gdt_2 + self.gdt_4) / 4
        return self._gdtha

    @property
    def rmsd(self):
        """ RMSD

        Computed on :attr:`~rigid_transformed_mapped_model_pos` and
        :attr:`~rigid_mapped_target_pos`

        :type: :class:`float`
        """
        if self._rmsd is None:
            LogScript("Computing RMSD")
            self._rmsd = \
            self.rigid_mapped_target_pos.GetRMSD(self.rigid_transformed_mapped_model_pos)
        return self._rmsd

    @property
    def cad_score(self):
        """ The global CAD atom-atom (AA) score

        Computed based on :attr:`~model`. In case of oligomers, :attr:`~mapping`
        is used.

        :type: :class:`float`
        """
        if self._cad_score is None:
            self._compute_cad_score()
        return self._cad_score

    @property
    def local_cad_score(self):
        """ The per-residue CAD atom-atom (AA) scores

        Computed based on :attr:`~model`. In case of oligomers, :attr:`~mapping`
        is used.

        :type: :class:`dict`
        """
        if self._local_cad_score is None:
            self._compute_cad_score()
        return self._local_cad_score

    @property
    def patch_qs(self):
        """ Patch QS-scores for each residue in :attr:`~model_interface_residues`

        Representative patches for each residue r in chain c are computed as
        follows:
    
        * mdl_patch_one: All residues in c with CB (CA for GLY) positions within
          8A of r and within 12A of residues from any other chain.
        * mdl_patch_two: Closest residue x to r in any other chain gets
          identified. Patch is then constructed by selecting all residues from
          any other chain within 8A of x and within 12A from any residue in c.
        * trg_patch_one: Chain name and residue number based mapping from
          mdl_patch_one
        * trg_patch_two: Chain name and residue number based mapping from
          mdl_patch_two

        Results are stored in the same manner as
        :attr:`~model_interface_residues`, with corresponding scores instead of
        residue numbers. Scores for residues which are not
        :class:`mol.ChemType.AMINOACIDS` are set to None. Additionally,
        interface patches are derived from :attr:`~model`. If they contain
        residues which are not covered by :attr:`~target`, the score is set to
        None too.

        :type: :class:`dict` with chain names as key and and :class:`list`
                with scores of the respective interface residues.
        """
        if self._patch_qs is None:
            self._compute_patchqs_scores()
        return self._patch_qs

    @property
    def patch_dockq(self):
        """ Same as :attr:`~patch_qs` but for DockQ scores
        """
        if self._patch_dockq is None:
            self._compute_patchdockq_scores()
        return self._patch_dockq

    @property
    def tm_score(self):
        """ TM-score computed with USalign

        USalign executable can be specified with usalign_exec kwarg at Scorer
        construction, an OpenStructure internal copy of the USalign code is
        used otherwise.

        :type: :class:`float`
        """
        if self._tm_score is None:
            self._compute_tmscore()
        return self._tm_score

    @property
    def usalign_mapping(self):
        """ Mapping computed with USalign

        Dictionary with target chain names as key and model chain names as
        values. No guarantee that all chains are mapped. USalign executable
        can be specified with usalign_exec kwarg at Scorer construction, an
        OpenStructure internal copy of the USalign code is used otherwise.

        :type: :class:`dict`
        """
        if self._usalign_mapping is None:
            self._compute_tmscore()
        return self._usalign_mapping

    def _aln_helper(self, target, model):
        # perform required alignments - cannot take the alignments from the
        # mapping results as we potentially remove stuff there as compared
        # to self.model and self.target
        trg_seqs = dict()
        for ch in target.chains:
            cname = ch.GetName()
            s = ''.join([r.one_letter_code for r in ch.residues])
            s = seq.CreateSequence(ch.GetName(), s)
            trg_seqs[ch.GetName()] = s
        mdl_seqs = dict()
        for ch in model.chains:
            cname = ch.GetName()
            s = ''.join([r.one_letter_code for r in ch.residues])
            s = seq.CreateSequence(cname, s)
            mdl_seqs[ch.GetName()] = s

        alns = list()
        trg_pep_chains = [s.GetName() for s in self.chain_mapper.polypep_seqs]
        trg_nuc_chains = [s.GetName() for s in self.chain_mapper.polynuc_seqs]
        trg_pep_chains = set(trg_pep_chains)
        trg_nuc_chains = set(trg_nuc_chains)
        for trg_ch, mdl_ch in self.mapping.GetFlatMapping().items():
            if mdl_ch in mdl_seqs and trg_ch in trg_seqs:
                if trg_ch in trg_pep_chains:
                    stype = mol.ChemType.AMINOACIDS
                elif trg_ch in trg_nuc_chains:
                    stype = mol.ChemType.NUCLEOTIDES
                else:
                    raise RuntimeError("Chain name inconsistency... ask "
                                       "Gabriel")
                if self.resnum_alignments:
                    aln = self.chain_mapper.ResNumAlign(trg_seqs[trg_ch],
                                                        mdl_seqs[mdl_ch],
                                                        target, model)
                else:
                    aln = self.chain_mapper.NWAlign(trg_seqs[trg_ch],
                                                    mdl_seqs[mdl_ch],
                                                    stype)

                alns.append(aln)
        return alns

    def _compute_aln(self):
        self._aln = self._aln_helper(self.target, self.model)

    def _compute_stereochecked_aln(self):
        # lets not redo the alignment and derive it from self.aln
        alns = list()
        for a in self.aln:
            trg_s = a.GetSequence(0)
            mdl_s = a.GetSequence(1)
            trg_ch = self.target.FindChain(trg_s.name)
            mdl_ch = self.model.FindChain(mdl_s.name)

            sc_trg_olc = ['-'] * len(trg_s)
            sc_mdl_olc = ['-'] * len(mdl_s)

            sc_trg_ch = self.stereochecked_target.FindChain(trg_s.name)
            if sc_trg_ch.IsValid():
                # there is the theoretical possibility that the full chain
                # has been removed in stereochemistry checks...
                trg_residues = trg_ch.residues
                res_idx = 0
                for olc_idx, olc in enumerate(trg_s):
                    if olc != '-':
                        r = trg_residues[res_idx]
                        sc_r = sc_trg_ch.FindResidue(r.GetNumber())
                        if sc_r.IsValid():
                            sc_trg_olc[olc_idx] = sc_r.one_letter_code
                        res_idx += 1

            sc_mdl_ch = self.stereochecked_model.FindChain(mdl_s.name)
            if sc_mdl_ch.IsValid():
                # there is the theoretical possibility that the full chain
                # has been removed in stereochemistry checks...
                mdl_residues = mdl_ch.residues
                res_idx = 0
                for olc_idx, olc in enumerate(mdl_s):
                    if olc != '-':
                        r = mdl_residues[res_idx]
                        sc_r = sc_mdl_ch.FindResidue(r.GetNumber())
                        if sc_r.IsValid():
                            sc_mdl_olc[olc_idx] = sc_r.one_letter_code
                        res_idx += 1

            sc_trg_s = seq.CreateSequence(trg_s.name, ''.join(sc_trg_olc))
            sc_mdl_s = seq.CreateSequence(mdl_s.name, ''.join(sc_mdl_olc))
            new_a = seq.CreateAlignment()
            new_a.AddSequence(sc_trg_s)
            new_a.AddSequence(sc_mdl_s)
            alns.append(new_a)

        self._stereochecked_aln = alns

    def _compute_pepnuc_aln(self):
        self._pepnuc_aln = self._aln_helper(self.pepnuc_target,
                                            self.pepnuc_model)

    def _compute_lddt(self):
        LogScript("Computing all-atom LDDT")
        # LDDT requires a flat mapping with mdl_ch as key and trg_ch as value
        flat_mapping = self.mapping.GetFlatMapping(mdl_as_key=True)

        # make alignments accessible by mdl seq name
        alns = dict()
        for aln in self.aln:
            mdl_seq = aln.GetSequence(1)
            alns[mdl_seq.name] = aln

        # score variables to be set
        lddt_score = None
        local_lddt = None
        aa_local_lddt = None

        if self.lddt_no_stereochecks:
            lddt_chain_mapping = dict()
            for mdl_ch, trg_ch in flat_mapping.items():
                if mdl_ch in alns:
                    lddt_chain_mapping[mdl_ch] = trg_ch
            lddt_score = self.lddt_scorer.lDDT(self.model,
                                               chain_mapping = lddt_chain_mapping,
                                               residue_mapping = alns,
                                               check_resnames=False,
                                               local_lddt_prop="lddt",
                                               add_mdl_contacts = self.lddt_add_mdl_contacts,
                                               set_atom_props=True)[0]
            local_lddt = dict()
            aa_local_lddt = dict()
            for r in self.model.residues:

                cname = r.GetChain().GetName()
                if cname not in local_lddt:
                    local_lddt[cname] = dict()
                    aa_local_lddt[cname] = dict()

                rnum = r.GetNumber()
                if rnum not in aa_local_lddt[cname]:
                    aa_local_lddt[cname][rnum] = dict()

                if r.HasProp("lddt"):
                    score = round(r.GetFloatProp("lddt"), 3)
                    local_lddt[cname][rnum] = score
                else:
                    # not covered by trg or skipped in chain mapping procedure
                    # the latter happens if its part of a super short chain
                    local_lddt[cname][rnum] = None

                for a in r.atoms:
                    if a.HasProp("lddt"):
                        score = round(a.GetFloatProp("lddt"), 3)
                        aa_local_lddt[cname][rnum][a.GetName()] = score
                    else:
                        # not covered by trg or skipped in chain mapping
                        # procedure the latter happens if its part of a
                        # super short chain
                        aa_local_lddt[cname][rnum][a.GetName()] = None


        else:
            # keep track what chains we have in the stereochecked model
            # there might be really wild cases where a full model
            # chain is removed in the stereochemistry checks.
            # We need to adapt lddt chain mapping in these cases
            mdl_chains = set([ch.name for ch in self.stereochecked_model.chains])

            # make alignments accessible by mdl seq name
            stereochecked_alns = dict()
            for aln in self.stereochecked_aln:
                mdl_seq = aln.GetSequence(1)
                if mdl_seq.GetName() in mdl_chains:
                    stereochecked_alns[mdl_seq.name] = aln

            lddt_chain_mapping = dict()
            for mdl_ch, trg_ch in flat_mapping.items():
                if mdl_ch in stereochecked_alns:
                    lddt_chain_mapping[mdl_ch] = trg_ch


            lddt_score = self.lddt_scorer.lDDT(self.stereochecked_model,
                                               chain_mapping = lddt_chain_mapping,
                                               residue_mapping = stereochecked_alns,
                                               check_resnames=False,
                                               local_lddt_prop="lddt",
                                               add_mdl_contacts = self.lddt_add_mdl_contacts,
                                               set_atom_props=True)[0]
            local_lddt = dict()
            aa_local_lddt = dict()
            for r in self.model.residues:

                cname = r.GetChain().GetName()
                if cname not in local_lddt:
                    local_lddt[cname] = dict()
                    aa_local_lddt[cname] = dict()

                rnum = r.GetNumber()
                if rnum not in aa_local_lddt[cname]:
                    aa_local_lddt[cname][rnum] = dict()

                if r.HasProp("lddt"):
                    score = round(r.GetFloatProp("lddt"), 3)
                    local_lddt[cname][rnum] = score

                    trg_r = None # represents stereochecked trg res
                    mdl_r = None # represents stereochecked mdl res

                    for a in r.atoms:
                        if a.HasProp("lddt"):
                            score = round(a.GetFloatProp("lddt"), 3)
                            aa_local_lddt[cname][rnum][a.GetName()] = score
                        else:
                            # the target residue is there since we have a score
                            # for the residue.
                            # opt 1: The atom was never there in the
                            #        stereochecked target => None
                            # opt 2: The atom has been removed in the model
                            #        stereochecks but is there in stereochecked
                            #        target => 0.0
                            if trg_r is None:
                                # let's first see if we find that target residue
                                # in the non-stereochecked target
                                tmp = None
                                if cname in flat_mapping:
                                    for x, y in _GetAlignedResidues(alns[cname],
                                                                    self.target,
                                                                    self.model):
                                        if y.number == r.number:
                                            tmp = x
                                            break
                                if tmp is not None:
                                    # we have it in the non-stereochecked target!
                                    tmp_cname = tmp.GetChain().GetName()
                                    tmp_rnum = tmp.GetNumber()
                                    tmp = self.stereochecked_target.FindResidue(tmp_cname,
                                                                                tmp_rnum)
                                    if tmp.IsValid():
                                        # And it's there in the stereochecked target too!
                                        trg_r = tmp

                            if mdl_r is None:
                                tmp = self.stereochecked_model.FindResidue(cname, rnum)
                                if tmp.IsValid():
                                    mdl_r = tmp

                            if trg_r is None:
                                # opt 1 - the whole target residue is not there
                                # this is actually an impossibility, as we have
                                # a score for the full mdl residue set
                                aa_local_lddt[cname][rnum][a.GetName()] = None                                
                            elif trg_r is not None and not trg_r.FindAtom(a.GetName()).IsValid():
                                # opt 1 - the target residue is there but not the atom
                                aa_local_lddt[cname][rnum][a.GetName()] = None
                            elif trg_r is not None and trg_r.FindAtom(a.GetName()).IsValid() and \
                                 mdl_r is None:
                                # opt 2 - trg atom is there but full model residue is removed
                                # this is actuall an impossibility, as we have
                                # a score for the full mdl residue set
                                aa_local_lddt[cname][rnum][a.GetName()] = 0.0
                            elif trg_r is not None and trg_r.FindAtom(a.GetName()).IsValid() and \
                                 mdl_r is not None and not mdl_r.FindAtom(a.GetName()).IsValid():
                                # opt 2 - trg atom is there but model atom is removed
                                aa_local_lddt[cname][rnum][a.GetName()] = 0.0
                            else:
                                # unknown issue
                                aa_local_lddt[cname][rnum][a.GetName()] = None

                else:
                    mdl_res = self.stereochecked_model.FindResidue(cname, rnum)
                    if mdl_res.IsValid():
                        # not covered by trg or skipped in chain mapping procedure
                        # the latter happens if its part of a super short chain
                        local_lddt[cname][rnum] = None
                        for a in r.atoms:
                            aa_local_lddt[cname][rnum][a.GetName()] = None
                    else:
                        # opt 1: removed by stereochecks => assign 0.0
                        # opt 2: removed by stereochecks AND not covered by ref
                        #        => assign None

                        # fetch trg residue from non-stereochecked aln
                        trg_r = None
                        if cname in flat_mapping:
                            tmp = None
                            for x, y in _GetAlignedResidues(alns[cname],
                                                            self.target,
                                                            self.model):
                                if y.number == r.number:
                                    tmp = x
                                    break
                            if tmp is not None:
                                # we have it in the non-stereochecked target!
                                tmp_cname = tmp.GetChain().GetName()
                                tmp_rnum = tmp.GetNumber()
                                tmp = self.stereochecked_target.FindResidue(tmp_cname,
                                                                            tmp_rnum)
                                if tmp.IsValid():
                                    # And it's there in the stereochecked target too!
                                    trg_r = tmp

                        if trg_r is None:
                            local_lddt[cname][rnum] = None
                            for a in r.atoms:
                                aa_local_lddt[cname][rnum][a.GetName()] = None
                        else:
                            local_lddt[cname][rnum] = 0.0
                            for a in r.atoms:
                                if trg_r.FindAtom(a.GetName()).IsValid():
                                    aa_local_lddt[cname][rnum][a.GetName()] = 0.0
                                else:
                                    aa_local_lddt[cname][rnum][a.GetName()] = None

        self._lddt = lddt_score
        self._local_lddt = local_lddt
        self._aa_local_lddt = aa_local_lddt

    def _compute_bb_lddt(self):
        LogScript("Computing backbone LDDT")
        # make alignments accessible by mdl seq name
        alns = dict()
        for aln in self.aln:
            mdl_seq = aln.GetSequence(1)
            alns[mdl_seq.name] = aln

        # LDDT requires a flat mapping with mdl_ch as key and trg_ch as value
        flat_mapping = self.mapping.GetFlatMapping(mdl_as_key=True)
        lddt_chain_mapping = dict()
        for mdl_ch, trg_ch in flat_mapping.items():
            if mdl_ch in alns:
                lddt_chain_mapping[mdl_ch] = trg_ch

        lddt_score = self.bb_lddt_scorer.lDDT(self.model,
                                              chain_mapping = lddt_chain_mapping,
                                              residue_mapping = alns,
                                              check_resnames=False,
                                              local_lddt_prop="bb_lddt",
                                              add_mdl_contacts = self.lddt_add_mdl_contacts)[0]
        local_lddt = dict()
        for r in self.model.residues:
            cname = r.GetChain().GetName()
            if cname not in local_lddt:
                local_lddt[cname] = dict()
            if r.HasProp("bb_lddt"):
                score = round(r.GetFloatProp("bb_lddt"), 3)
                local_lddt[cname][r.GetNumber()] = score
            else:
                # not covered by trg or skipped in chain mapping procedure
                # the latter happens if its part of a super short chain
                local_lddt[cname][r.GetNumber()] = None

        self._bb_lddt = lddt_score
        self._bb_local_lddt = local_lddt

    def _compute_ilddt(self):
        LogScript("Computing all-atom iLDDT")
        # LDDT requires a flat mapping with mdl_ch as key and trg_ch as value
        flat_mapping = self.mapping.GetFlatMapping(mdl_as_key=True)

        if self.lddt_no_stereochecks:
            alns = dict()
            for aln in self.aln:
                mdl_seq = aln.GetSequence(1)
                alns[mdl_seq.name] = aln
            lddt_chain_mapping = dict()
            for mdl_ch, trg_ch in flat_mapping.items():
                if mdl_ch in alns:
                    lddt_chain_mapping[mdl_ch] = trg_ch
            self._ilddt = self.lddt_scorer.lDDT(self.model,
                                                chain_mapping = lddt_chain_mapping,
                                                residue_mapping = alns,
                                                check_resnames=False,
                                                local_lddt_prop="lddt",
                                                add_mdl_contacts = self.lddt_add_mdl_contacts,
                                                no_intrachain=True)[0]
        else:

            # keep track what chains we have in the stereochecked model
            # there might be really wild cases where a full model
            # chain is removed in the stereochemistry checks.
            # We need to adapt lddt chain mapping in these cases
            mdl_chains = set([ch.name for ch in self.stereochecked_model.chains])

            # make alignments accessible by mdl seq name
            stereochecked_alns = dict()
            for aln in self.stereochecked_aln:
                mdl_seq = aln.GetSequence(1)
                if mdl_seq.GetName() in mdl_chains:
                    stereochecked_alns[mdl_seq.name] = aln

            lddt_chain_mapping = dict()
            for mdl_ch, trg_ch in flat_mapping.items():
                if mdl_ch in stereochecked_alns:
                    lddt_chain_mapping[mdl_ch] = trg_ch

            self._ilddt = self.lddt_scorer.lDDT(self.stereochecked_model,
                                                chain_mapping = lddt_chain_mapping,
                                                residue_mapping = stereochecked_alns,
                                                check_resnames=False,
                                                local_lddt_prop="lddt",
                                                add_mdl_contacts = self.lddt_add_mdl_contacts,
                                                no_intrachain=True)[0]


    def _compute_qs(self):
        LogScript("Computing global QS-score")
        qs_score_result = self.qs_scorer.Score(self.mapping.mapping)
        self._qs_global = qs_score_result.QS_global
        self._qs_best = qs_score_result.QS_best

    def _compute_per_interface_qs_scores(self):
        LogScript("Computing per-interface QS-score")
        self._per_interface_qs_global = list()
        self._per_interface_qs_best = list()

        for interface in self.qs_interfaces:
            trg_ch1 = interface[0]
            trg_ch2 = interface[1]
            mdl_ch1 = interface[2]
            mdl_ch2 = interface[3]
            qs_res = self.qs_scorer.ScoreInterface(trg_ch1, trg_ch2,
                                                   mdl_ch1, mdl_ch2)
            self._per_interface_qs_best.append(qs_res.QS_best)
            self._per_interface_qs_global.append(qs_res.QS_global)

    def _compute_ics_scores(self):
        LogScript("Computing ICS scores")
        contact_scorer_res = self.contact_scorer.ScoreICS(self.mapping.mapping)
        self._ics_precision = contact_scorer_res.precision
        self._ics_recall = contact_scorer_res.recall
        self._ics = contact_scorer_res.ics
        self._per_interface_ics_precision = list()
        self._per_interface_ics_recall = list()
        self._per_interface_ics = list()
        flat_mapping = self.mapping.GetFlatMapping()
        for trg_int in self.contact_target_interfaces:
            trg_ch1 = trg_int[0]
            trg_ch2 = trg_int[1]
            if trg_ch1 in flat_mapping and trg_ch2 in flat_mapping:
                mdl_ch1 = flat_mapping[trg_ch1]
                mdl_ch2 = flat_mapping[trg_ch2]
                res = self.contact_scorer.ScoreICSInterface(trg_ch1, trg_ch2,
                                                            mdl_ch1, mdl_ch2)
                self._per_interface_ics_precision.append(res.precision)
                self._per_interface_ics_recall.append(res.recall)
                self._per_interface_ics.append(res.ics)
            else:
                self._per_interface_ics_precision.append(None)
                self._per_interface_ics_recall.append(None)
                self._per_interface_ics.append(None)

    def _trim_model(self):
        trimmed_mdl = mol.CreateEntityFromView(self.mapping.model, True)
        trimmed_aln = dict()

        for trg_cname, mdl_cname in self.mapping.GetFlatMapping().items():
            mdl_ch = trimmed_mdl.FindChain(mdl_cname)
            aln = self.mapping.alns[(trg_cname, mdl_cname)]

            # some limited test that stuff matches
            assert(mdl_ch.GetResidueCount() == \
                   len(aln.GetSequence(1).GetGaplessString()))

            mdl_residues = mdl_ch.residues
            mdl_res_idx = 0
            aligned_mdl_seq = ['-'] * aln.GetLength()
            for col_idx, col in enumerate(aln):
                if col[0] != '-' and col[1] != '-':
                    mdl_res = mdl_residues[mdl_res_idx]
                    mdl_res.SetIntProp("aligned", 1)
                    aligned_mdl_seq[col_idx] = col[1]
                if col[1] != '-':
                    mdl_res_idx += 1
            aligned_mdl_seq = ''.join(aligned_mdl_seq)
            aligned_mdl_seq = seq.CreateSequence(mdl_cname, aligned_mdl_seq)

            new_aln = seq.CreateAlignment()
            new_aln.AddSequence(aln.GetSequence(0))
            new_aln.AddSequence(aligned_mdl_seq)
            trimmed_aln[(trg_cname, mdl_cname)] = new_aln

        self._trimmed_model = trimmed_mdl.Select("graligned:0=1")
        self._trimmed_aln = trimmed_aln

    def _compute_ics_scores_trimmed(self):
        LogScript("Computing ICS scores trimmed")

        # this is an ugly hack without any efficiency in mind
        # we're simply taking the entities from mapper and construct
        # a new contact scorer from scratch

        contact_scorer_res = self.trimmed_contact_scorer.ScoreICS(self.mapping.mapping)
        self._ics_trimmed = contact_scorer_res.ics
        self._ics_precision_trimmed = contact_scorer_res.precision
        self._ics_recall_trimmed = contact_scorer_res.recall

        self._per_interface_ics_precision_trimmed = list()
        self._per_interface_ics_recall_trimmed = list()
        self._per_interface_ics_trimmed = list()
        flat_mapping = self.mapping.GetFlatMapping()
        for trg_int in self.contact_target_interfaces:
            trg_ch1 = trg_int[0]
            trg_ch2 = trg_int[1]
            if trg_ch1 in flat_mapping and trg_ch2 in flat_mapping:
                mdl_ch1 = flat_mapping[trg_ch1]
                mdl_ch2 = flat_mapping[trg_ch2]
                res = self.trimmed_contact_scorer.ScoreICSInterface(trg_ch1, trg_ch2,
                                                                    mdl_ch1, mdl_ch2)
                self._per_interface_ics_precision_trimmed.append(res.precision)
                self._per_interface_ics_recall_trimmed.append(res.recall)
                self._per_interface_ics_trimmed.append(res.ics)
            else:
                self._per_interface_ics_precision_trimmed.append(None)
                self._per_interface_ics_recall_trimmed.append(None)
                self._per_interface_ics_trimmed.append(None)

    def _compute_ips_scores(self):
        LogScript("Computing IPS scores")
        contact_scorer_res = self.contact_scorer.ScoreIPS(self.mapping.mapping)
        self._ips_precision = contact_scorer_res.precision
        self._ips_recall = contact_scorer_res.recall
        self._ips = contact_scorer_res.ips

        self._per_interface_ips_precision = list()
        self._per_interface_ips_recall = list()
        self._per_interface_ips = list()
        flat_mapping = self.mapping.GetFlatMapping()
        for trg_int in self.contact_target_interfaces:
            trg_ch1 = trg_int[0]
            trg_ch2 = trg_int[1]
            if trg_ch1 in flat_mapping and trg_ch2 in flat_mapping:
                mdl_ch1 = flat_mapping[trg_ch1]
                mdl_ch2 = flat_mapping[trg_ch2]
                res = self.contact_scorer.ScoreIPSInterface(trg_ch1, trg_ch2,
                                                            mdl_ch1, mdl_ch2)
                self._per_interface_ips_precision.append(res.precision)
                self._per_interface_ips_recall.append(res.recall)
                self._per_interface_ips.append(res.ips)
            else:
                self._per_interface_ips_precision.append(None)
                self._per_interface_ips_recall.append(None)
                self._per_interface_ips.append(None)

    def _compute_ips_scores_trimmed(self):
        LogScript("Computing IPS scores trimmed")

        # this is an ugly hack without any efficiency in mind
        # we're simply taking the entities from mapper and construct
        # a new contact scorer from scratch
        contact_scorer_res = self.trimmed_contact_scorer.ScoreIPS(self.mapping.mapping)
        self._ips_precision_trimmed = contact_scorer_res.precision
        self._ips_recall_trimmed = contact_scorer_res.recall
        self._ips_trimmed = contact_scorer_res.ips

        self._per_interface_ips_precision_trimmed = list()
        self._per_interface_ips_recall_trimmed = list()
        self._per_interface_ips_trimmed = list()
        flat_mapping = self.mapping.GetFlatMapping()
        for trg_int in self.contact_target_interfaces:
            trg_ch1 = trg_int[0]
            trg_ch2 = trg_int[1]
            if trg_ch1 in flat_mapping and trg_ch2 in flat_mapping:
                mdl_ch1 = flat_mapping[trg_ch1]
                mdl_ch2 = flat_mapping[trg_ch2]
                res = self.trimmed_contact_scorer.ScoreIPSInterface(trg_ch1, trg_ch2,
                                                                    mdl_ch1, mdl_ch2)
                self._per_interface_ips_precision_trimmed.append(res.precision)
                self._per_interface_ips_recall_trimmed.append(res.recall)
                self._per_interface_ips_trimmed.append(res.ips)
            else:
                self._per_interface_ips_precision_trimmed.append(None)
                self._per_interface_ips_recall_trimmed.append(None)
                self._per_interface_ips_trimmed.append(None)

    def _compute_dockq_scores(self):
        LogScript("Computing DockQ")

        if self.dockq_capri_peptide and len(self.chain_mapper.polynuc_seqs) > 0:
            raise RuntimeError("Cannot compute DockQ for reference structures "
                               "with nucleotide chains if dockq_capri_peptide "
                               "is enabled.")

        # lists with values in contact_target_interfaces
        self._dockq_scores = list()
        self._fnat = list()
        self._fnonnat = list()
        self._irmsd = list()
        self._lrmsd = list()
        self._nnat = list()
        self._nmdl = list()

        dockq_alns = dict()
        for aln in self.aln:
            trg_s = aln.GetSequence(0)
            mdl_s = aln.GetSequence(1)
            dockq_alns[(trg_s.GetName(), mdl_s.GetName())] = aln

        for interface in self.dockq_interfaces:
            trg_ch1 = interface[0]
            trg_ch2 = interface[1]
            mdl_ch1 = interface[2]
            mdl_ch2 = interface[3]
            aln1 = dockq_alns[(trg_ch1, mdl_ch1)]
            aln2 = dockq_alns[(trg_ch2, mdl_ch2)]
            if self.dockq_capri_peptide:
                res = dockq.DockQ(self.model, self.target, mdl_ch1, mdl_ch2,
                                  trg_ch1, trg_ch2, ch1_aln=aln1,
                                  ch2_aln=aln2, contact_dist_thresh = 4.0,
                                  interface_dist_thresh=8.0,
                                  interface_cb = True)
            else:
                res = dockq.DockQ(self.model, self.target, mdl_ch1, mdl_ch2,
                                  trg_ch1, trg_ch2, ch1_aln=aln1,
                                  ch2_aln=aln2)

            self._fnat.append(res["fnat"])
            self._fnonnat.append(res["fnonnat"])
            self._irmsd.append(res["irmsd"])
            self._lrmsd.append(res["lrmsd"])
            self._dockq_scores.append(res["DockQ"])
            self._nnat.append(res["nnat"])
            self._nmdl.append(res["nmdl"])

        # keep track of native counts in target interfaces which are
        # not covered in model in order to compute
        # dockq_ave_full/dockq_wave_full in the end
        not_covered_counts = list()
        proc_trg_interfaces = set([(x[0], x[1]) for x in self.dockq_interfaces])
        for interface in self.dockq_target_interfaces:
            if interface not in proc_trg_interfaces:
                # let's run DockQ with trg as trg/mdl in order to get the native
                # contacts out - no need to pass alns as the residue numbers
                # match for sure
                trg_ch1 = interface[0]
                trg_ch2 = interface[1]

                if self.dockq_capri_peptide:
                    res = dockq.DockQ(self.target, self.target,
                                      trg_ch1, trg_ch2, trg_ch1, trg_ch2,
                                      contact_dist_thresh = 4.0,
                                      interface_dist_thresh=8.0,
                                      interface_cb = True)
                else:
                    res = dockq.DockQ(self.target, self.target,
                                      trg_ch1, trg_ch2, trg_ch1, trg_ch2)

                not_covered_counts.append(res["nnat"])
  
        # there are 4 types of combined scores
        # - simple average
        # - average weighted by native_contacts
        # - the two above including nonmapped_contact_interfaces => set DockQ to 0.0
        scores = np.array(self._dockq_scores)
        weights = np.array(self._nnat)
        if len(scores) > 0:
            self._dockq_ave = np.mean(scores)
        else:
            self._dockq_ave = 0.0
        self._dockq_wave = np.sum(np.multiply(weights/np.sum(weights), scores))
        scores = np.append(scores, [0.0]*len(not_covered_counts))
        weights = np.append(weights, not_covered_counts)
        if len(scores) > 0:
            self._dockq_ave_full = np.mean(scores)
        else:
            self._dockq_ave_full = 0.0
        self._dockq_wave_full = np.sum(np.multiply(weights/np.sum(weights),
                                                   scores))

    def _extract_mapped_pos(self):
        self._mapped_target_pos = geom.Vec3List()
        self._mapped_model_pos = geom.Vec3List()
        self._n_target_not_mapped = 0
        processed_trg_chains = set()
        for trg_ch, mdl_ch in self.mapping.GetFlatMapping().items():
            processed_trg_chains.add(trg_ch)
            aln = self.mapping.alns[(trg_ch, mdl_ch)]
            n_mapped = 0
            for trg_res, mdl_res in _GetAlignedResidues(aln,
                                                        self.mapping.target,
                                                        self.mapping.model):
                trg_at = trg_res.FindAtom("CA")
                mdl_at = mdl_res.FindAtom("CA")
                if not trg_at.IsValid():
                    trg_at = trg_res.FindAtom("C3'")
                if not mdl_at.IsValid():
                    mdl_at = mdl_res.FindAtom("C3'")
                self._mapped_target_pos.append(trg_at.GetPos())
                self._mapped_model_pos.append(mdl_at.GetPos())
                n_mapped += 1
            self._n_target_not_mapped += (len(aln.GetSequence(0).GetGaplessString())-n_mapped)
        # count number of trg residues from non-mapped chains
        for ch in self.mapping.target.chains:
            if ch.GetName() not in processed_trg_chains:
                self._n_target_not_mapped += ch.GetResidueCount()

    def _extract_mapped_pos_full_bb(self):
        self._mapped_target_pos_full_bb = geom.Vec3List()
        self._mapped_model_pos_full_bb = geom.Vec3List()
        exp_pep_atoms = ["N", "CA", "C"]
        exp_nuc_atoms = ["O5'", "C5'", "C4'", "C3'", "O3'"]
        trg_pep_chains = [s.GetName() for s in self.chain_mapper.polypep_seqs]
        trg_nuc_chains = [s.GetName() for s in self.chain_mapper.polynuc_seqs]
        for trg_ch, mdl_ch in self.mapping.GetFlatMapping().items():
            aln = self.mapping.alns[(trg_ch, mdl_ch)]
            trg_ch = aln.GetSequence(0).GetName()
            if trg_ch in trg_pep_chains:
                exp_atoms = exp_pep_atoms
            elif trg_ch in trg_nuc_chains:
                exp_atoms = exp_nuc_atoms
            else:
                # this should be guaranteed by the chain mapper
                raise RuntimeError("Unexpected error - contact OST developer")
            for trg_res, mdl_res in _GetAlignedResidues(aln,
                                                        self.mapping.target,
                                                        self.mapping.model):
                for aname in exp_atoms:
                    trg_at = trg_res.FindAtom(aname)
                    mdl_at = mdl_res.FindAtom(aname)
                    if not (trg_at.IsValid() and mdl_at.IsValid()):
                        # this should be guaranteed by the chain mapper
                        raise RuntimeError("Unexpected error - contact OST "
                                           "developer")
                    self._mapped_target_pos_full_bb.append(trg_at.GetPos())
                    self._mapped_model_pos_full_bb.append(mdl_at.GetPos())


    def _extract_rigid_mapped_pos(self):
        self._rigid_mapped_target_pos = geom.Vec3List()
        self._rigid_mapped_model_pos = geom.Vec3List()
        self._rigid_n_target_not_mapped = 0
        processed_trg_chains = set()
        for trg_ch, mdl_ch in self.rigid_mapping.GetFlatMapping().items():
            processed_trg_chains.add(trg_ch)
            aln = self.rigid_mapping.alns[(trg_ch, mdl_ch)]
            n_mapped = 0
            for trg_res, mdl_res in _GetAlignedResidues(aln,
                                                        self.rigid_mapping.target,
                                                        self.rigid_mapping.model):
                trg_at = trg_res.FindAtom("CA")
                mdl_at = mdl_res.FindAtom("CA")
                if not trg_at.IsValid():
                    trg_at = trg_res.FindAtom("C3'")
                if not mdl_at.IsValid():
                    mdl_at = mdl_res.FindAtom("C3'")
                self._rigid_mapped_target_pos.append(trg_at.GetPos())
                self._rigid_mapped_model_pos.append(mdl_at.GetPos())
                n_mapped += 1

            self._rigid_n_target_not_mapped += (len(aln.GetSequence(0).GetGaplessString())-n_mapped)
        # count number of trg residues from non-mapped chains
        for ch in self.rigid_mapping.target.chains:
            if ch.GetName() not in processed_trg_chains:
                self._rigid_n_target_not_mapped += ch.GetResidueCount()

    def _extract_rigid_mapped_pos_full_bb(self):
        self._rigid_mapped_target_pos_full_bb = geom.Vec3List()
        self._rigid_mapped_model_pos_full_bb = geom.Vec3List()
        exp_pep_atoms = ["N", "CA", "C"]
        exp_nuc_atoms = ["O5'", "C5'", "C4'", "C3'", "O3'"]
        trg_pep_chains = [s.GetName() for s in self.chain_mapper.polypep_seqs]
        trg_nuc_chains = [s.GetName() for s in self.chain_mapper.polynuc_seqs]
        for trg_ch, mdl_ch in self.rigid_mapping.GetFlatMapping().items():
            aln = self.rigid_mapping.alns[(trg_ch, mdl_ch)]
            trg_ch = aln.GetSequence(0).GetName()
            if trg_ch in trg_pep_chains:
                exp_atoms = exp_pep_atoms
            elif trg_ch in trg_nuc_chains:
                exp_atoms = exp_nuc_atoms
            else:
                # this should be guaranteed by the chain mapper
                raise RuntimeError("Unexpected error - contact OST developer")
            for trg_res, mdl_res in _GetAlignedResidues(aln,
                                                        self.rigid_mapping.target,
                                                        self.rigid_mapping.model):
                for aname in exp_atoms:
                    trg_at = trg_res.FindAtom(aname)
                    mdl_at = mdl_res.FindAtom(aname)
                    if not (trg_at.IsValid() and mdl_at.IsValid()):
                        # this should be guaranteed by the chain mapper
                        raise RuntimeError("Unexpected error - contact OST "
                                           "developer")
                    self._rigid_mapped_target_pos_full_bb.append(trg_at.GetPos())
                    self._rigid_mapped_model_pos_full_bb.append(mdl_at.GetPos())

    def _compute_cad_score(self):
        if not self.resnum_alignments:
            raise RuntimeError("CAD score computations rely on residue numbers "
                               "that are consistent between target and model "
                               "chains, i.e. only work if resnum_alignments "
                               "is True at Scorer construction.")
        try:
            LogScript("Computing CAD score")
            cad_score_exec = \
            settings.Locate("voronota-cadscore",
                            explicit_file_name=self.cad_score_exec)
        except Exception as e:
            raise RuntimeError("voronota-cadscore must be in PATH for CAD "
                               "score scoring") from e
        cad_bin_dir = os.path.dirname(cad_score_exec)
        m = self.mapping.GetFlatMapping(mdl_as_key=True)
        cad_result = cadscore.CADScore(self.model, self.target,
                                       mode = "voronota",
                                       label="localcad",
                                       old_regime=False,
                                       cad_bin_path=cad_bin_dir,
                                       chain_mapping=m)

        local_cad = dict()
        for r in self.model.residues:
            cname = r.GetChain().GetName()
            if cname not in local_cad:
                local_cad[cname] = dict()
            if r.HasProp("localcad"):
                score = round(r.GetFloatProp("localcad"), 3)
                local_cad[cname][r.GetNumber()] = score
            else:
                local_cad[cname][r.GetNumber()] = None

        self._cad_score = cad_result.globalAA
        self._local_cad_score = local_cad

    def _get_repr_view(self, ent):
        """ Returns view with representative peptide atoms => CB, CA for GLY
    
        Ensures that each residue has exactly one atom with assertions
    
        :param ent: Entity for which you want the representative view
        :param ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        :returns: :class:`ost.mol.EntityView` derived from ent
        """
        repr_ent = ent.Select("(aname=\"CB\" or (rname=\"GLY\" and aname=\"CA\"))")
        for r in repr_ent.residues:
            assert(len(r.atoms) == 1)
        return repr_ent

    def _get_interface_residues(self, ent):
        """ Get interface residues
    
        Thats all residues having a contact with at least one residue from another
        chain (CB-CB distance <= 8A, CA in case of Glycine)
    
        :param ent: Model for which to extract interface residues
        :type ent: :class:`ost.mol.EntityView`
        :returns: :class:`dict` with chain names as key and and :class:`list`
                  with residue numbers of the respective interface residues.
        """
        # select for representative positions => CB, CA for GLY 
        repr_ent = self._get_repr_view(ent)
        result = {ch.GetName(): list() for ch in ent.chains}
        for ch in ent.chains:
            cname = ch.GetName()
            sel = repr_ent.Select(f"(cname={mol.QueryQuoteName(cname)} and 8 <> [cname!={mol.QueryQuoteName(cname)}])")
            result[cname] = [r.GetNumber() for r in sel.residues]
        return result

    def _do_stereochecks(self):
        """ Perform stereochemistry checks on model and target
        """
        LogInfo("Performing stereochemistry checks on model and target")
        data = stereochemistry.GetDefaultStereoData()
        l_data = stereochemistry.GetDefaultStereoLinkData()

        a, b, c, d = stereochemistry.StereoCheck(self.model, stereo_data = data,
                                                 stereo_link_data = l_data)
        self._stereochecked_model = a
        self._model_clashes = b
        self._model_bad_bonds = c
        self._model_bad_angles = d

        a, b, c, d = stereochemistry.StereoCheck(self.target, stereo_data = data,
                                                 stereo_link_data = l_data)

        self._stereochecked_target = a
        self._target_clashes = b
        self._target_bad_bonds = c
        self._target_bad_angles = d

    def _get_interface_patches(self, mdl_ch, mdl_rnum):
        """ Select interface patches representative for specified residue
    
        The patches for specified residue r in chain c are selected as follows:
    
        * mdl_patch_one: All residues in c with CB (CA for GLY) positions within 8A
          of r and within 12A of residues from any other chain.
        * mdl_patch_two: Closest residue x to r in any other chain gets identified.
          Patch is then constructed by selecting all residues from any other chain
          within 8A of x and within 12A from any residue in c.
        * trg_patch_one: Chain name and residue number based mapping from
          mdl_patch_one
        * trg_patch_two: Chain name and residue number based mapping from
          mdl_patch_two
    
        :param mdl_ch: Name of chain in *self.model* of residue of interest
        :type mdl_ch: :class:`str`
        :param mdl_rnum: Residue number of residue of interest
        :type mdl_rnum: :class:`ost.mol.ResNum`
        :returns: Tuple with 5 elements: 1) :class:`bool` flag whether all residues
                  in *mdl* patches are covered in *trg* 2) mtl_patch_one
                  3) mdl_patch_two 4) trg_patch_one 5) trg_patch_two
        """
        # select for representative positions => CB, CA for GLY 
        repr_mdl = self._get_repr_view(self.model.Select("peptide=true"))
    
        # get position for specified residue
        r = self.model.FindResidue(mdl_ch, mdl_rnum)
        if not r.IsValid():
            raise RuntimeError(f"Cannot find residue {mdl_rnum} in chain {mdl_ch}")
        if r.GetName() == "GLY":
            at = r.FindAtom("CA")
        else:
            at = r.FindAtom("CB")
        if not at.IsValid():
            raise RuntimeError("Cannot find interface views for res without CB/CA")
        r_pos = at.GetPos()
    
        # mdl_patch_one contains residues from the same chain as r
        # => all residues within 8A of r and within 12A of any other chain
    
        # q1 selects for everything in same chain and within 8A of r_pos
        q1 = f"(cname={mol.QueryQuoteName(mdl_ch)} and 8 <> {{{r_pos[0]},{r_pos[1]},{r_pos[2]}}})"
        # q2 selects for everything within 12A of any other chain
        q2 = f"(12 <> [cname!={mol.QueryQuoteName(mdl_ch)}])"
        mdl_patch_one = self.model.CreateEmptyView()
        sel = repr_mdl.Select(" and ".join([q1, q2]))
        for r in sel.residues:
            mdl_r = self.model.FindResidue(r.GetChain().GetName(), r.GetNumber())
            mdl_patch_one.AddResidue(mdl_r, mol.ViewAddFlag.INCLUDE_ALL)
    
        # mdl_patch_two contains residues from all other chains. In detail:
        # the closest residue to r is identified in any other chain, and the
        # patch is filled with residues that are within 8A of that residue and
        # within 12A of chain from r
        sel = repr_mdl.Select(f"(cname!={mol.QueryQuoteName(mdl_ch)})")
        close_stuff = sel.FindWithin(r_pos, 8)
        min_pos = None
        min_dist = 42.0
        for close_at in close_stuff:
            dist = geom.Distance(r_pos, close_at.GetPos())
            if dist < min_dist:
                min_pos = close_at.GetPos()
                min_dist = dist
    
        # q1 selects for everything not in mdl_ch but within 8A of min_pos
        q1 = f"(cname!={mol.QueryQuoteName(mdl_ch)} and 8 <> {{{min_pos[0]},{min_pos[1]},{min_pos[2]}}})"
        # q2 selects for everything within 12A of mdl_ch
        q2 = f"(12 <> [cname={mol.QueryQuoteName(mdl_ch)}])"
        mdl_patch_two = self.model.CreateEmptyView()
        sel = repr_mdl.Select(" and ".join([q1, q2]))
        for r in sel.residues:
            mdl_r = self.model.FindResidue(r.GetChain().GetName(), r.GetNumber())
            mdl_patch_two.AddResidue(mdl_r, mol.ViewAddFlag.INCLUDE_ALL)
    
        # transfer mdl residues to trg
        flat_mapping = self.mapping.GetFlatMapping(mdl_as_key=True)
        full_trg_coverage = True
        trg_patch_one = self.mapping.target.CreateEmptyView()
        for r in mdl_patch_one.residues:
            trg_r = None
            mdl_cname = r.GetChain().GetName()
            if mdl_cname in flat_mapping:
                aln = self.mapping.alns[(flat_mapping[mdl_cname], mdl_cname)]
                for x,y in _GetAlignedResidues(aln,
                                               self.mapping.target,
                                               self.mapping.model):
                    if y.GetNumber() == r.GetNumber():
                        trg_r = x
                        break

            if trg_r is not None:
                trg_patch_one.AddResidue(trg_r.handle,
                                         mol.ViewAddFlag.INCLUDE_ALL)
            else:
                full_trg_coverage = False
    
        trg_patch_two = self.mapping.target.CreateEmptyView()
        for r in mdl_patch_two.residues:
            trg_r = None
            mdl_cname = r.GetChain().GetName()
            if mdl_cname in flat_mapping:
                aln = self.mapping.alns[(flat_mapping[mdl_cname], mdl_cname)]
                for x,y in _GetAlignedResidues(aln,
                                               self.mapping.target,
                                               self.mapping.model):
                    if y.GetNumber() == r.GetNumber():
                        trg_r = x
                        break

            if trg_r is not None:
                trg_patch_two.AddResidue(trg_r.handle,
                                         mol.ViewAddFlag.INCLUDE_ALL)
            else:
                full_trg_coverage = False
    
        return (full_trg_coverage, mdl_patch_one, mdl_patch_two,
                trg_patch_one, trg_patch_two)

    def _compute_patchqs_scores(self):
        LogScript("Computing patch QS-scores")
        self._patch_qs = dict()
        for cname, rnums in self.model_interface_residues.items():
            scores = list()
            for rnum in rnums:
                score = None
                r = self.model.FindResidue(cname, rnum)
                if r.IsValid() and r.GetChemType() == mol.ChemType.AMINOACIDS:
                    full_trg_coverage, mdl_patch_one, mdl_patch_two, \
                    trg_patch_one, trg_patch_two = \
                    self._get_interface_patches(cname, rnum)
                    if full_trg_coverage:
                        score = self._patchqs(mdl_patch_one, mdl_patch_two,
                                              trg_patch_one, trg_patch_two)
                scores.append(score)
            self._patch_qs[cname] = scores

    def _compute_patchdockq_scores(self):
        LogScript("Computing patch DockQ scores")
        self._patch_dockq = dict()
        for cname, rnums in self.model_interface_residues.items():
            scores = list()
            for rnum in rnums:
                score = None
                r = self.model.FindResidue(cname, rnum)
                if r.IsValid() and r.GetChemType() == mol.ChemType.AMINOACIDS:
                    full_trg_coverage, mdl_patch_one, mdl_patch_two, \
                    trg_patch_one, trg_patch_two = \
                    self._get_interface_patches(cname, rnum)
                    if full_trg_coverage:
                        score = self._patchdockq(mdl_patch_one, mdl_patch_two,
                                                 trg_patch_one, trg_patch_two)
                scores.append(score)
            self._patch_dockq[cname] = scores

    def _patchqs(self, mdl_patch_one, mdl_patch_two, trg_patch_one, trg_patch_two):
        """ Score interface residue patches with QS-score
    
        In detail: Construct two entities with two chains each. First chain
        consists of residues from <x>_patch_one and second chain consists of
        <x>_patch_two. The returned score is the QS-score between the two
        entities
    
        :param mdl_patch_one: Interface patch representing scored residue
        :type mdl_patch_one: :class:`ost.mol.EntityView`
        :param mdl_patch_two: Interface patch representing scored residue
        :type mdl_patch_two: :class:`ost.mol.EntityView`
        :param trg_patch_one: Interface patch representing scored residue
        :type trg_patch_one: :class:`ost.mol.EntityView`
        :param trg_patch_two: Interface patch representing scored residue
        :type trg_patch_two: :class:`ost.mol.EntityView`
        :returns: PatchQS score
        """
        qs_ent_mdl = self._qs_ent_from_patches(mdl_patch_one, mdl_patch_two)
        qs_ent_trg = self._qs_ent_from_patches(trg_patch_one, trg_patch_two)
    
        alnA = seq.CreateAlignment()
        s = ''.join([r.one_letter_code for r in mdl_patch_one.residues])
        alnA.AddSequence(seq.CreateSequence("A", s))
        s = ''.join([r.one_letter_code for r in trg_patch_one.residues])
        alnA.AddSequence(seq.CreateSequence("A", s))
    
        alnB = seq.CreateAlignment()
        s = ''.join([r.one_letter_code for r in mdl_patch_two.residues])
        alnB.AddSequence(seq.CreateSequence("B", s))
        s = ''.join([r.one_letter_code for r in trg_patch_two.residues])
        alnB.AddSequence(seq.CreateSequence("B", s))
        alns = {("A", "A"): alnA, ("B", "B"): alnB}
    
        scorer = QSScorer(qs_ent_mdl, [["A"], ["B"]], qs_ent_trg, alns)
        score_result = scorer.Score([["A"], ["B"]])
    
        return score_result.QS_global
    
    def _patchdockq(self, mdl_patch_one, mdl_patch_two, trg_patch_one,
                    trg_patch_two):
        """ Score interface residue patches with DockQ
    
        In detail: Construct two entities with two chains each. First chain
        consists of residues from <x>_patch_one and second chain consists of
        <x>_patch_two. The returned score is the QS-score between the two
        entities
    
        :param mdl_patch_one: Interface patch representing scored residue
        :type mdl_patch_one: :class:`ost.mol.EntityView`
        :param mdl_patch_two: Interface patch representing scored residue
        :type mdl_patch_two: :class:`ost.mol.EntityView`
        :param trg_patch_one: Interface patch representing scored residue
        :type trg_patch_one: :class:`ost.mol.EntityView`
        :param trg_patch_two: Interface patch representing scored residue
        :type trg_patch_two: :class:`ost.mol.EntityView`
        :returns: DockQ score
        """
        m = self._qs_ent_from_patches(mdl_patch_one, mdl_patch_two)
        t = self._qs_ent_from_patches(trg_patch_one, trg_patch_two)
        dockq_result = dockq.DockQ(t, m, "A", "B", "A", "B")
        if dockq_result["nnat"] > 0:
            return dockq_result["DockQ"]
        return 0.0

    def _qs_ent_from_patches(self, patch_one, patch_two):
        """ Constructs Entity with two chains named "A" and "B""
    
        Blindly adds all residues from *patch_one* to chain A and residues from
        patch_two to chain B.
        """
        ent = mol.CreateEntity()
        ed = ent.EditXCS()
        added_ch = ed.InsertChain("A")
        for r in patch_one.residues:
            added_r = ed.AppendResidue(added_ch, r.GetName())
            added_r.SetChemClass(str(r.GetChemClass()))
            for a in r.atoms:
                ed.InsertAtom(added_r, a.handle)
        added_ch = ed.InsertChain("B")
        for r in patch_two.residues:
            added_r = ed.AppendResidue(added_ch, r.GetName())
            added_r.SetChemClass(str(r.GetChemClass()))
            for a in r.atoms:
                ed.InsertAtom(added_r, a.handle)
        return ent

    def _construct_custom_mapping(self, mapping):
        """ constructs a full blown MappingResult object from simple dict

        :param mapping: mapping with trg chains as key and mdl ch as values
        :type mapping: :class:`dict`
        """
        LogInfo("Setting custom chain mapping")

        chain_mapper = self.chain_mapper
        chem_mapping, chem_group_alns, mdl_chains_without_chem_mapping, mdl = \
        chain_mapper.GetChemMapping(self.model)

        # now that we have a chem mapping, lets do consistency checks
        # - check whether chain names are unique and available in structures
        # - check whether the mapped chains actually map to the same chem groups
        if len(mapping) != len(set(mapping.keys())):
            raise RuntimeError(f"Expected unique trg chain names in mapping. Got "
                               f"{mapping.keys()}")
        if len(mapping) != len(set(mapping.values())):
            raise RuntimeError(f"Expected unique mdl chain names in mapping. Got "
                               f"{mapping.values()}")

        trg_chains = set([ch.GetName() for ch in chain_mapper.target.chains])
        mdl_chains = set([ch.GetName() for ch in mdl.chains])
        for k,v in mapping.items():
            if k not in trg_chains:
                raise RuntimeError(f"Target chain \"{k}\" is not available "
                                   f"in target processed for chain mapping "
                                   f"({trg_chains})")
            if v not in mdl_chains:
                raise RuntimeError(f"Model chain \"{v}\" is not available "
                                   f"in model processed for chain mapping "
                                   f"({mdl_chains})")

        for trg_ch, mdl_ch in mapping.items():
            trg_group_idx = None
            mdl_group_idx = None
            for idx, group in enumerate(chain_mapper.chem_groups):
                if trg_ch in group:
                    trg_group_idx = idx
                    break
            for idx, group in enumerate(chem_mapping):
                if mdl_ch in group:
                    mdl_group_idx = idx
                    break
            if trg_group_idx is None or mdl_group_idx is None:
                raise RuntimeError("Could not establish a valid chem grouping "
                                   "of chain names provided in custom mapping.")
            
            if trg_group_idx != mdl_group_idx:
                raise RuntimeError(f"Chem group mismatch in custom mapping: "
                                   f"target chain \"{trg_ch}\" groups with the "
                                   f"following chemically equivalent target "
                                   f"chains: "
                                   f"{chain_mapper.chem_groups[trg_group_idx]} "
                                   f"but model chain \"{mdl_ch}\" maps to the "
                                   f"following target chains: "
                                   f"{chain_mapper.chem_groups[mdl_group_idx]}")

        pairs = set([(trg_ch, mdl_ch) for trg_ch, mdl_ch in mapping.items()])
        ref_mdl_alns =  \
        chain_mapping._GetRefMdlAlns(chain_mapper.chem_groups,
                                     chain_mapper.chem_group_alignments,
                                     chem_mapping,
                                     chem_group_alns,
                                     pairs = pairs)

        # translate mapping format
        final_mapping = list()
        for ref_chains in chain_mapper.chem_groups:
            mapped_mdl_chains = list()
            for ref_ch in ref_chains:
                if ref_ch in mapping:
                    mapped_mdl_chains.append(mapping[ref_ch])
                else:
                    mapped_mdl_chains.append(None)
            final_mapping.append(mapped_mdl_chains)

        alns = dict()
        for ref_group, mdl_group in zip(chain_mapper.chem_groups,
                                        final_mapping):
            for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                if ref_ch is not None and mdl_ch is not None:
                    aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                    trg_view = chain_mapper.target.Select(f"cname={mol.QueryQuoteName(ref_ch)}")
                    mdl_view = mdl.Select(f"cname={mol.QueryQuoteName(mdl_ch)}")
                    aln.AttachView(0, trg_view)
                    aln.AttachView(1, mdl_view)
                    alns[(ref_ch, mdl_ch)] = aln

        return chain_mapping.MappingResult(chain_mapper.target, mdl,
                                           chain_mapper.chem_groups,
                                           chem_mapping,
                                           mdl_chains_without_chem_mapping,
                                           final_mapping, alns)

    def _compute_tmscore(self):
        res = None
        if self.usalign_exec is None:
            LogScript("Computing TM-score with built-in USalign")
            if self.oum:
                flat_mapping = self.rigid_mapping.GetFlatMapping()
                LogInfo("Overriding TM-score chain mapping")
                res = res = bindings.WrappedMMAlign(self.model, self.target,
                                                    mapping=flat_mapping)
            else:
                res = bindings.WrappedMMAlign(self.model, self.target)
        else:
            LogScript("Computing TM-score with USalign exectuable")
            if self.oum:
                LogInfo("Overriding TM-score chain mapping")
                flat_mapping = self.rigid_mapping.GetFlatMapping()
                res = tmtools.USAlign(self.model, self.target,
                                      usalign = self.usalign_exec,
                                      custom_chain_mapping = flat_mapping)
            else:
                res = tmtools.USAlign(self.model, self.target,
                                      usalign = self.usalign_exec)

        self._tm_score = res.tm_score
        self._usalign_mapping = dict()
        for a,b in zip(res.ent1_mapped_chains, res.ent2_mapped_chains):
            self._usalign_mapping[b] = a

    def _resnum_connect(self, ent):
        ed = None
        for ch in ent.chains:
            res_list = ch.residues
            for i in range(len(res_list) - 1):
                ra = res_list[i]
                rb = res_list[i+1]
                if ra.GetNumber().GetNum() + 1 == rb.GetNumber().GetNum():
                    if ra.IsPeptideLinking() and rb.IsPeptideLinking():
                        c = ra.FindAtom("C")
                        n = rb.FindAtom("N")
                        if c.IsValid() and n.IsValid() and not mol.BondExists(c, n):
                            if ed is None:
                                ed = ent.EditXCS(mol.BUFFERED_EDIT)
                            ed.Connect(c,n,1)
                    elif ra.IsNucleotideLinking() and rb.IsNucleotideLinking():
                        o = ra.FindAtom("O3'")
                        p = rb.FindAtom("P")
                        if o.IsValid() and p.IsValid()and not mol.BondExists(o, p):
                            if ed is None:
                                ed = ent.EditXCS(mol.BUFFERED_EDIT)
                            ed.Connect(o,p,1)


# specify public interface
__all__ = ('lDDTBSScorer', 'Scorer',)
