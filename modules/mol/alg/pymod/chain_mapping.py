"""
Chain mapping aims to identify a one-to-one relationship between chains in a
reference structure and a model.
"""

import itertools
import copy

import numpy as np

from scipy.special import factorial
from scipy.special import binom # as of Python 3.8, the math module implements
                                # comb, i.e. n choose k

import ost
from ost import seq
from ost import mol
from ost import geom

from ost.mol.alg import lddt
from ost.mol.alg import bb_lddt
from ost.mol.alg import qsscore

def _CSel(ent, cnames):
    """ Returns view with specified chains

    Ensures that quotation marks are around chain names to not confuse
    OST query language with weird special characters.
    """
    query = "cname=" + ','.join([mol.QueryQuoteName(cname) for cname in cnames])
    return ent.Select(query)

class MappingResult:
    """ Result object for the chain mapping functions in :class:`ChainMapper`

    Constructor is directly called within the functions, no need to construct
    such objects yourself.
    """
    def __init__(self, target, model, chem_groups, chem_mapping,
                 unmapped_mdl_chains, mapping, alns, opt_score=None):
        self._target = target
        self._model = model
        self._chem_groups = chem_groups
        self._chem_mapping = chem_mapping
        self._unmapped_mdl_chains = unmapped_mdl_chains
        self._mapping = mapping
        self._alns = alns
        self._opt_score = opt_score

    @property
    def target(self):
        """ Target/reference structure, i.e. :attr:`ChainMapper.target`

        :type: :class:`ost.mol.EntityView`
        """
        return self._target

    @property
    def model(self):
        """ Model structure that gets mapped onto :attr:`~target`

        Underwent same processing as :attr:`ChainMapper.target`, i.e.
        only contains peptide/nucleotide chains of sufficient size.

        :type: :class:`ost.mol.EntityView`
        """
        return self._model

    @property
    def chem_groups(self):
        """ Groups of chemically equivalent chains in :attr:`~target`

        Same as :attr:`ChainMapper.chem_group`

        :class:`list` of :class:`list` of :class:`str` (chain names)
        """
        return self._chem_groups

    @property
    def chem_mapping(self):
        """ Assigns chains in :attr:`~model` to :attr:`~chem_groups`.

        :class:`list` of :class:`list` of :class:`str` (chain names)
        """
        return self._chem_mapping

    @property
    def unmapped_mdl_chains(self):
        """ Model chains that cannot be mapped to :attr:`chem_groups`

        Depends on parameterization of :class:`ChainMapper`

        :class:`list` of class:`str` (chain names)
        """
        return self._unmapped_mdl_chains

    @property
    def mapping(self):
        """ Mapping of :attr:`~model` chains onto :attr:`~target`

        Exact same shape as :attr:`~chem_groups` but containing the names of the
        mapped chains in :attr:`~model`. May contain None for :attr:`~target`
        chains that are not covered. No guarantee that all chains in
        :attr:`~model` are mapped.

        :class:`list` of :class:`list` of :class:`str` (chain names)
        """
        return self._mapping

    @property
    def alns(self):
        """ Alignments of mapped chains in :attr:`~target` and :attr:`~model`

        Each alignment is accessible with ``alns[(t_chain,m_chain)]``. First
        sequence is the sequence of :attr:`target` chain, second sequence the
        one from :attr:`~model`.

        :type: :class:`dict` with key: :class:`tuple` of :class:`str`, value:
               :class:`ost.seq.AlignmentHandle`
        """
        return self._alns

    @property
    def opt_score(self):
        """ Placeholder property without any guarantee of being set

        Different scores get optimized in the various chain mapping algorithms.
        Some of them may set their final optimal score in that property.
        Consult the documentation of the respective chain mapping algorithm
        for more information. Won't be in the return dict of
        :func:`JSONSummary`.
        """
        return self._opt_score

    def GetFlatMapping(self, mdl_as_key=False):
        """ Returns flat mapping as :class:`dict` for all mapable chains

        :param mdl_as_key: Default is target chain name as key and model chain
                           name as value. This can be reversed with this flag.
        :returns: :class:`dict` with :class:`str` as key/value that describe
                  one-to-one mapping
        """
        flat_mapping = dict()
        for trg_chem_group, mdl_chem_group in zip(self.chem_groups,
                                                  self.mapping):
            for a,b in zip(trg_chem_group, mdl_chem_group):
                if a is not None and b is not None:
                    if mdl_as_key:
                        flat_mapping[b] = a
                    else:
                        flat_mapping[a] = b
        return flat_mapping

    def JSONSummary(self):
        """ Returns JSON serializable summary of results
        """
        json_dict = dict()
        json_dict["chem_groups"] = self.chem_groups
        json_dict["mapping"] = self.mapping
        json_dict["flat_mapping"] = self.GetFlatMapping()
        json_dict["alns"] = list()
        for aln in self.alns.values():
            trg_seq = aln.GetSequence(0)
            mdl_seq = aln.GetSequence(1)
            aln_dict = {"trg_ch": trg_seq.GetName(), "trg_seq": str(trg_seq),
                        "mdl_ch": mdl_seq.GetName(), "mdl_seq": str(mdl_seq)}
            json_dict["alns"].append(aln_dict)
        return json_dict


class ReprResult:

    """ Result object for :func:`ChainMapper.GetRepr`

    Constructor is directly called within the function, no need to construct
    such objects yourself.

    :param lDDT: lDDT for this mapping. Depends on how you call
                 :func:`ChainMapper.GetRepr` whether this is backbone only or
                 full atom lDDT.
    :type lDDT: :class:`float`
    :param substructure: The full substructure for which we searched for a
                         representation
    :type substructure: :class:`ost.mol.EntityView`
    :param ref_view: View pointing to the same underlying entity as
                     *substructure* but only contains the stuff that is mapped
    :type ref_view: :class:`mol.EntityView`
    :param mdl_view: The matching counterpart in model
    :type mdl_view: :class:`mol.EntityView`
    """
    def __init__(self, lDDT, substructure, ref_view, mdl_view):
        self._lDDT = lDDT
        self._substructure = substructure
        assert(len(ref_view.residues) == len(mdl_view.residues))
        self._ref_view = ref_view
        self._mdl_view = mdl_view

        # lazily evaluated attributes
        self._ref_bb_pos = None
        self._mdl_bb_pos = None
        self._ref_full_bb_pos = None
        self._mdl_full_bb_pos = None
        self._superposition = None
        self._superposed_mdl_bb_pos = None
        self._ost_query = None
        self._flat_mapping = None
        self._inconsistent_residues = None

    @property
    def lDDT(self):
        """ lDDT of representation result

        Depends on how you call :func:`ChainMapper.GetRepr` whether this is
        backbone only or full atom lDDT.

        :type: :class:`float`
        """
        return self._lDDT

    @property
    def substructure(self):
        """ The full substructure for which we searched for a
        representation

        :type: :class:`ost.mol.EntityView`
        """
        return self._substructure

    @property
    def ref_view(self):
        """ View which contains the mapped subset of :attr:`substructure`

        :type: :class:`ost.mol.EntityView`
        """
        return self._ref_view

    @property
    def mdl_view(self):
        """ The :attr:`ref_view` representation in the model

        :type: :class:`ost.mol.EntityView`
        """
        return self._mdl_view
    
    @property
    def ref_residues(self):
        """ The reference residues

        :type: class:`mol.ResidueViewList`
        """
        return self.ref_view.residues
    
    @property
    def mdl_residues(self):
        """ The model residues

        :type: :class:`mol.ResidueViewList`
        """
        return self.mdl_view.residues

    @property
    def inconsistent_residues(self):
        """ A list of mapped residue whose names do not match (eg. ALA in the
        reference and LEU in the model).

        The mismatches are reported as a tuple of :class:`~ost.mol.ResidueView`
        (reference, model), or as an empty list if all the residue names match.

        :type: :class:`list`
        """
        if self._inconsistent_residues is None:
            self._inconsistent_residues = self._GetInconsistentResidues(
                self.ref_residues, self.mdl_residues)
        return self._inconsistent_residues

    @property
    def ref_bb_pos(self):
        """ Representative backbone positions for reference residues.

        Thats CA positions for peptides and C3' positions for Nucleotides.

        :type: :class:`geom.Vec3List`
        """
        if self._ref_bb_pos is None:
            self._ref_bb_pos = self._GetBBPos(self.ref_residues)
        return self._ref_bb_pos

    @property
    def mdl_bb_pos(self):
        """ Representative backbone positions for model residues.

        Thats CA positions for peptides and C3' positions for Nucleotides.

        :type: :class:`geom.Vec3List`
        """
        if self._mdl_bb_pos is None:
            self._mdl_bb_pos = self._GetBBPos(self.mdl_residues)
        return self._mdl_bb_pos

    @property
    def ref_full_bb_pos(self):
        """ Representative backbone positions for reference residues.

        Thats N, CA and C positions for peptides and O5', C5', C4', C3', O3'
        positions for Nucleotides.

        :type: :class:`geom.Vec3List`
        """
        if self._ref_full_bb_pos is None:
            self._ref_full_bb_pos = self._GetFullBBPos(self.ref_residues)
        return self._ref_full_bb_pos

    @property
    def mdl_full_bb_pos(self):
        """ Representative backbone positions for reference residues.

        Thats N, CA and C positions for peptides and O5', C5', C4', C3', O3'
        positions for Nucleotides.

        :type: :class:`geom.Vec3List`
        """
        if self._mdl_full_bb_pos is None:
            self._mdl_full_bb_pos = self._GetFullBBPos(self.mdl_residues)
        return self._mdl_full_bb_pos


    @property
    def superposition(self):
        """ Superposition of mdl residues onto ref residues

        Superposition computed as minimal RMSD superposition on
        :attr:`ref_bb_pos` and :attr:`mdl_bb_pos`. If number of positions is
        smaller 3, the full_bb_pos equivalents are used instead.

        :type: :class:`ost.mol.alg.SuperpositionResult`
        """
        if self._superposition is None:
            if len(self.mdl_bb_pos) < 3:
                self._superposition = _GetSuperposition(self.mdl_full_bb_pos,
                                                self.ref_full_bb_pos, False)
            else:
                self._superposition = _GetSuperposition(self.mdl_bb_pos,
                                                self.ref_bb_pos, False)
        return self._superposition

    @property
    def transform(self):
        """ Transformation to superpose mdl residues onto ref residues

        :type: :class:`ost.geom.Mat4`
        """
        return self.superposition.transformation

    @property
    def superposed_mdl_bb_pos(self):
        """ :attr:`mdl_bb_pos` with :attr:`transform applied`

        :type: :class:`geom.Vec3List`
        """
        if self._superposed_mdl_bb_pos is None:
            self._superposed_mdl_bb_pos = geom.Vec3List(self.mdl_bb_pos)
            self._superposed_mdl_bb_pos.ApplyTransform(self.transform)
        return self._superposed_mdl_bb_pos

    @property
    def bb_rmsd(self):
        """ RMSD of the binding site backbone atoms after :attr:`superposition`

        :type: :class:`float`
        """
        return self.superposition.rmsd


    @property
    def ost_query(self):
        """ query for mdl residues in OpenStructure query language

        Repr can be selected as ``full_mdl.Select(ost_query)``

        Returns invalid query if residue numbers have insertion codes.

        :type: :class:`str`
        """
        if self._ost_query is None:
            chain_rnums = dict()
            for r in self.mdl_residues:
                chname = r.GetChain().GetName()
                rnum = r.GetNumber().GetNum()
                if chname not in chain_rnums:
                    chain_rnums[chname] = list()
                chain_rnums[chname].append(str(rnum))
            chain_queries = list()
            for k,v in chain_rnums.items():
                q = f"(cname={mol.QueryQuoteName(k)} and "
                q += f"rnum={','.join(v)})"
                chain_queries.append(q)
            self._ost_query = " or ".join(chain_queries)
        return self._ost_query

    def JSONSummary(self):
        """ Returns JSON serializable summary of results
        """
        json_dict = dict()
        json_dict["lDDT"] = self.lDDT
        json_dict["ref_residues"] = [r.GetQualifiedName() for r in \
                                     self.ref_residues]
        json_dict["mdl_residues"] = [r.GetQualifiedName() for r in \
                                     self.mdl_residues]
        json_dict["transform"] = list(self.transform.data)
        json_dict["bb_rmsd"] = self.bb_rmsd
        json_dict["ost_query"] = self.ost_query
        json_dict["flat_mapping"] = self.GetFlatChainMapping()
        return json_dict

    def GetFlatChainMapping(self, mdl_as_key=False):
        """ Returns flat mapping of all chains in the representation

        :param mdl_as_key: Default is target chain name as key and model chain
                           name as value. This can be reversed with this flag.
        :returns: :class:`dict` with :class:`str` as key/value that describe
                  one-to-one mapping
        """
        flat_mapping = dict()
        for trg_res, mdl_res in zip(self.ref_residues, self.mdl_residues):
            if mdl_as_key:
                flat_mapping[mdl_res.chain.name] = trg_res.chain.name
            else:
                flat_mapping[trg_res.chain.name] = mdl_res.chain.name
        return flat_mapping

    def _GetFullBBPos(self, residues):
        """ Helper to extract full backbone positions
        """
        exp_pep_atoms = ["N", "CA", "C"]
        exp_nuc_atoms = ["O5'", "C5'", "C4'", "C3'", "O3'"]
        bb_pos = geom.Vec3List()
        for r in residues:
            if r.GetChemType() == mol.ChemType.NUCLEOTIDES:
                exp_atoms = exp_nuc_atoms
            elif r.GetChemType() == mol.ChemType.AMINOACIDS:
                exp_atoms = exp_pep_atoms
            else:
                raise RuntimeError(f"Residue {r.qualified_name} has unexpected"
                                   f" chem type {r.GetChemType()}")
            for aname in exp_atoms:
                a = r.FindAtom(aname)
                if not a.IsValid():
                    raise RuntimeError(f"Expected atom {aname} missing from "
                                      f"{r.qualified_name}")
                bb_pos.append(a.GetPos())
        return bb_pos

    def _GetBBPos(self, residues):
        """ Helper to extract single representative position for each residue
        """
        bb_pos = geom.Vec3List()
        for r in residues:
            at = r.FindAtom("CA")
            if not at.IsValid():
                at = r.FindAtom("C3'")
            if not at.IsValid():
                raise RuntimeError(f"Residue {r.qualified_name} missing "
                                   f"expected CA/C3' atom")
            bb_pos.append(at.GetPos())
        return bb_pos

    def _GetInconsistentResidues(self, ref_residues, mdl_residues):
        """ Helper to extract a list of inconsistent residues.
        """
        if len(ref_residues) != len(mdl_residues):
            raise ValueError("Something terrible happened... Reference and "
                             "model lengths differ... RUN...")
        inconsistent_residues = list()
        for ref_residue, mdl_residue in zip(ref_residues, mdl_residues):
            if ref_residue.name != mdl_residue.name:
                inconsistent_residues.append((ref_residue, mdl_residue))
        return inconsistent_residues


class ChainMapper:
    """ Class to compute chain mappings

    All algorithms are performed on processed structures which fulfill
    criteria as given in constructor arguments (*min_pep_length*,
    "min_nuc_length") and only contain residues which have all required backbone
    atoms. for peptide residues thats N, CA, C and CB (no CB for GLY), for
    nucleotide residues thats O5', C5', C4', C3' and O3'.

    Chain mapping is a three step process:

    * Group chemically identical chains in *target* using pairwise
      alignments that are either computed with Needleman-Wunsch (NW) or
      simply derived from residue numbers (*resnum_alignments* flag).
      In case of NW, *pep_subst_mat*, *pep_gap_open* and *pep_gap_ext*
      and their nucleotide equivalents are relevant. Two chains are
      considered identical if they fulfill the *pep_seqid_thr* and
      have at least *min_pep_length* aligned residues. Same logic
      is applied for nucleotides with respective thresholds.

    * Map chains in an input model to these groups. Parameters for generating
      alignments are the same as above. Model chains are mapped to the chem
      group with highest sequence identity. Optionally, you can set a minimum
      sequence identity threshold with *mdl_map_pep_seqid_thr*/
      *mdl_map_nuc_seqid_thr* to avoid nonsensical mappings. You can either
      get the group mapping with :func:`GetChemMapping` or directly call
      one of the full fletched one-to-one chain mapping functions which
      execute that step internally.

    * Obtain one-to-one mapping for chains in an input model and
      *target* with one of the available mapping functions. Just to get an
      idea of complexity. If *target* and *model* are octamers, there are
      ``8! = 40320`` possible chain mappings.

    :param target: Target structure onto which models are mapped.
                   Computations happen on a selection only containing
                   polypeptides and polynucleotides.
    :type target: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param resnum_alignments: Use residue numbers instead of
                              Needleman-Wunsch to compute pairwise
                              alignments. Relevant for :attr:`~chem_groups` 
                              and related attributes.
    :type resnum_alignments: :class:`bool`
    :param pep_seqid_thr: Threshold used to decide when two chains are
                          identical. 95 percent tolerates the few mutations
                          crystallographers like to do.
    :type pep_seqid_thr:  :class:`float`
    :param nuc_seqid_thr: Nucleotide equivalent for *pep_seqid_thr*
    :type nuc_seqid_thr:  :class:`float`
    :param pep_subst_mat: Substitution matrix to align peptide sequences,
                          irrelevant if *resnum_alignments* is True,
                          defaults to seq.alg.BLOSUM62
    :type pep_subst_mat: :class:`ost.seq.alg.SubstWeightMatrix`
    :param pep_gap_open: Gap open penalty to align peptide sequences,
                         irrelevant if *resnum_alignments* is True
    :type pep_gap_open: :class:`int`
    :param pep_gap_ext: Gap extension penalty to align peptide sequences,
                        irrelevant if *resnum_alignments* is True
    :type pep_gap_ext: :class:`int`
    :param nuc_subst_mat: Nucleotide equivalent for *pep_subst_mat*,
                          defaults to seq.alg.NUC44
    :type nuc_subst_mat: :class:`ost.seq.alg.SubstWeightMatrix`
    :param nuc_gap_open: Nucleotide equivalent for *pep_gap_open*
    :type nuc_gap_open: :class:`int`
    :param nuc_gap_ext: Nucleotide equivalent for *pep_gap_ext*
    :type nuc_gap_ext: :class:`int`
    :param min_pep_length: Minimal number of residues for a peptide chain to be
                           considered in target and in models.
    :type min_pep_length: :class:`int`
    :param min_nuc_length: Minimal number of residues for a nucleotide chain to be
                           considered in target and in models.
    :type min_nuc_length: :class:`int` 
    :param n_max_naive: Max possible chain mappings that are enumerated in
                        :func:`~GetNaivelDDTMapping` /
                        :func:`~GetDecomposerlDDTMapping`. A
                        :class:`RuntimeError` is raised in case of bigger
                        complexity.
    :type n_max_naive: :class:`int`
    :param mdl_map_pep_seqid_thr: To avoid non-sensical mapping of model chains
                                  to reference chem groups. If a value larger
                                  than 0.0 is provided, minimal criteria for
                                  assignment becomes a sequence identity of
                                  *mdl_map_pep_seqid_thr* and at least
                                  *min_pep_length* aligned residues. If set to
                                  zero, it simply assigns a model chain to the
                                  chem group with highest sequence identity.
    :type mdl_map_pep_seqid_thr: :class:`float`
    :param mdl_map_nuc_seqid_thr: Nucleotide equivalent of
                                  *mdl_map_pep_seqid_thr*
    :type mdl_map_nuc_seqid_thr: :class:`float`
    """
    def __init__(self, target, resnum_alignments=False,
                 pep_seqid_thr = 95., nuc_seqid_thr = 95., 
                 pep_subst_mat = seq.alg.BLOSUM62, pep_gap_open = -11,
                 pep_gap_ext = -1, nuc_subst_mat = seq.alg.NUC44,
                 nuc_gap_open = -4, nuc_gap_ext = -4,
                 min_pep_length = 6, min_nuc_length = 4,
                 n_max_naive = 1e8,
                 mdl_map_pep_seqid_thr = 0.,
                 mdl_map_nuc_seqid_thr = 0.):

        # attributes
        self.resnum_alignments = resnum_alignments
        self.pep_seqid_thr = pep_seqid_thr
        self.nuc_seqid_thr = nuc_seqid_thr
        self.min_pep_length = min_pep_length
        self.min_nuc_length = min_nuc_length
        self.n_max_naive = n_max_naive
        self.mdl_map_pep_seqid_thr = mdl_map_pep_seqid_thr
        self.mdl_map_nuc_seqid_thr = mdl_map_nuc_seqid_thr
        self.pep_subst_mat = pep_subst_mat
        self.pep_gap_open = pep_gap_open
        self.pep_gap_ext = pep_gap_ext
        self.nuc_subst_mat = nuc_subst_mat
        self.nuc_gap_open = nuc_gap_open
        self.nuc_gap_ext = nuc_gap_ext

        # lazy computed attributes
        self._chem_groups = None
        self._chem_group_alignments = None
        self._chem_group_ref_seqs = None
        self._chem_group_types = None

        # target structure preprocessing
        self._target, self._polypep_seqs, self._polynuc_seqs = \
        self.ProcessStructure(target)

    @property
    def target(self):
        """Target structure that only contains peptides/nucleotides

        Contains only residues that have the backbone representatives
        (CA for peptide and C3' for nucleotides) to avoid ATOMSEQ alignment
        inconsistencies when switching between all atom and backbone only
        representations.

        :type: :class:`ost.mol.EntityView`
        """
        return self._target

    @property
    def polypep_seqs(self):
        """Sequences of peptide chains in :attr:`~target`

        :type: :class:`ost.seq.SequenceList`
        """
        return self._polypep_seqs

    @property
    def polynuc_seqs(self):
        """Sequences of nucleotide chains in :attr:`~target`

        :type: :class:`ost.seq.SequenceList`
        """
        return self._polynuc_seqs
    
    @property
    def chem_groups(self):
        """Groups of chemically equivalent chains in :attr:`~target`

        First chain in group is the one with longest sequence.
      
        :getter: Computed on first use (cached)
        :type: :class:`list` of :class:`list` of :class:`str` (chain names)
        """
        if self._chem_groups is None:
            self._ComputeChemGroups()
        return self._chem_groups
    
    @property
    def chem_group_alignments(self):
        """MSA for each group in :attr:`~chem_groups`

        The first sequence is the reference sequence.
        The subsequent sequences represent the ATOMSEQ sequences in
        :attr:`~target` in same order as in :attr:`~chem_groups`.

        :getter: Computed on first use (cached)
        :type: :class:`ost.seq.AlignmentList`
        """
        if self._chem_group_alignments is None:
            self._ComputeChemGroups()
        return self._chem_group_alignments

    @property
    def chem_group_ref_seqs(self):
        """Reference sequence for each group in :attr:`~chem_groups`

        :getter: Computed on first use (cached)
        :type: :class:`ost.seq.SequenceList`
        """
        if self._chem_group_ref_seqs is None:
            self._ComputeChemGroups()
        return self._chem_group_ref_seqs

    @property
    def chem_group_types(self):
        """ChemType of each group in :attr:`~chem_groups`

        Specifying if groups are poly-peptides/nucleotides, i.e. 
        :class:`ost.mol.ChemType.AMINOACIDS` or
        :class:`ost.mol.ChemType.NUCLEOTIDES`
        
        :getter: Computed on first use (cached)
        :type: :class:`list` of :class:`ost.mol.ChemType`
        """
        if self._chem_group_types is None:
            self._ComputeChemGroups()
        return self._chem_group_types
        
    def GetChemMapping(self, model):
        """Maps sequences in *model* to chem_groups of target

        :param model: Model from which to extract sequences, a
                      selection that only includes peptides and nucleotides
                      is performed and returned along other results.
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :returns: Tuple with two lists of length `len(self.chem_groups)`, a list
                  and an :class:`ost.mol.EntityView` representing *model*:
                  1) Each element is a :class:`list` with mdl chain names that
                  map to the chem group at that position.
                  2) Each element is a :class:`ost.seq.AlignmentList` aligning
                  these mdl chain sequences to the chem group ref sequences.
                  3) The model chains that could not be mapped to the reference
                  4) A selection of *model* that only contains polypeptides and
                  polynucleotides whose ATOMSEQ exactly matches the sequence
                  info in the returned alignments.
        """
        mdl, mdl_pep_seqs, mdl_nuc_seqs = self.ProcessStructure(model)
        mapping = [list() for x in self.chem_groups]
        unmapped_mdl_chains = list()
        alns = [seq.AlignmentList() for x in self.chem_groups]

        for s in mdl_pep_seqs:
            idx, aln = self._MapSequence(self.chem_group_ref_seqs, 
                                         self.chem_group_types,
                                         s, mol.ChemType.AMINOACIDS, mdl,
                                         seq_id_thr = self.mdl_map_pep_seqid_thr,
                                         min_aln_length = self.min_pep_length)
            if idx is None:
                unmapped_mdl_chains.append(s.GetName())
            else:
                mapping[idx].append(s.GetName())
                alns[idx].append(aln)

        for s in mdl_nuc_seqs:
            idx, aln = self._MapSequence(self.chem_group_ref_seqs, 
                                         self.chem_group_types,
                                         s, mol.ChemType.NUCLEOTIDES, mdl,
                                         seq_id_thr = self.mdl_map_nuc_seqid_thr,
                                         min_aln_length = self.min_nuc_length)
            if idx is None:
                unmapped_mdl_chains.append(s.GetName())
            else:
                mapping[idx].append(s.GetName())
                alns[idx].append(aln)

        return (mapping, alns, unmapped_mdl_chains, mdl)


    def GetlDDTMapping(self, model, inclusion_radius=15.0,
                       thresholds=[0.5, 1.0, 2.0, 4.0], strategy="heuristic",
                       steep_opt_rate = None, block_seed_size = 5,
                       block_blocks_per_chem_group = 5,
                       chem_mapping_result = None,
                       heuristic_n_max_naive = 40320):
        """ Identify chain mapping by optimizing lDDT score

        Maps *model* chain sequences to :attr:`~chem_groups` and find mapping
        based on backbone only lDDT score (CA for amino acids C3' for
        Nucleotides).

        Either performs a naive search, i.e. enumerate all possible mappings or
        executes a greedy strategy that tries to identify a (close to) optimal
        mapping in an iterative way by starting from a start mapping (seed). In
        each iteration, the one-to-one mapping that leads to highest increase
        in number of conserved contacts is added with the additional requirement
        that this added mapping must have non-zero interface counts towards the
        already mapped chains. So basically we're "growing" the mapped structure
        by only adding connected stuff.

        The available strategies:

        * **naive**: Enumerates all possible mappings and returns best        

        * **greedy_fast**: perform all vs. all single chain lDDTs within the
          respective ref/mdl chem groups. The mapping with highest number of
          conserved contacts is selected as seed for greedy extension

        * **greedy_full**: try multiple seeds for greedy extension, i.e. try
          all ref/mdl chain combinations within the respective chem groups and
          retain the mapping leading to the best lDDT.

        * **greedy_block**: try multiple seeds for greedy extension, i.e. try
          all ref/mdl chain combinations within the respective chem groups and
          extend them to *block_seed_size*. *block_blocks_per_chem_group*
          for each chem group are selected for exhaustive extension.

        * **heuristic**: Uses *naive* strategy if number of possible mappings
          is within *heuristic_n_max_naive*. The default of 40320 corresponds
          to an octamer (8!=40320). A structure with stoichiometry A6B2 would be
          6!*2!=1440 etc. If the number of possible mappings is larger,
          *greedy_full* is used.

        Sets :attr:`MappingResult.opt_score` in case of no trivial one-to-one
        mapping. 

        :param model: Model to map
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param inclusion_radius: Inclusion radius for lDDT
        :type inclusion_radius: :class:`float`
        :param thresholds: Thresholds for lDDT
        :type thresholds: :class:`list` of :class:`float`
        :param strategy: Strategy to find mapping. Must be in ["naive",
                         "greedy_fast", "greedy_full", "greedy_block"]
        :type strategy: :class:`str`
        :param steep_opt_rate: Only relevant for greedy strategies.
                               If set, every *steep_opt_rate* mappings, a simple
                               optimization is executed with the goal of
                               avoiding local minima. The optimization
                               iteratively checks all possible swaps of mappings
                               within their respective chem groups and accepts
                               swaps that improve lDDT score. Iteration stops as
                               soon as no improvement can be achieved anymore.
        :type steep_opt_rate: :class:`int`
        :param block_seed_size: Param for *greedy_block* strategy - Initial seeds
                                are extended by that number of chains.
        :type block_seed_size: :class:`int`
        :param block_blocks_per_chem_group: Param for *greedy_block* strategy -
                                            Number of blocks per chem group that
                                            are extended in an initial search
                                            for high scoring local solutions.
        :type block_blocks_per_chem_group: :class:`int`
        :param chem_mapping_result: Pro param. The result of
                                    :func:`~GetChemMapping` where you provided
                                    *model*. If set, *model* parameter is not
                                    used.
        :type chem_mapping_result: :class:`tuple`
        :returns: A :class:`MappingResult`
        """

        strategies = ["naive", "greedy_fast", "greedy_full", "greedy_block",
                      "heuristic"]
        if strategy not in strategies:
            raise RuntimeError(f"Strategy must be in {strategies}")

        if chem_mapping_result is None:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            self.GetChemMapping(model)
        else:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            chem_mapping_result

        ref_mdl_alns =  _GetRefMdlAlns(self.chem_groups,
                                       self.chem_group_alignments,
                                       chem_mapping,
                                       chem_group_alns)

        # check for the simplest case
        one_to_one = _CheckOneToOneMapping(self.chem_groups, chem_mapping)
        if one_to_one is not None:
            alns = dict()
            for ref_group, mdl_group in zip(self.chem_groups, one_to_one):
                for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                    if ref_ch is not None and mdl_ch is not None:
                        aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                        alns[(ref_ch, mdl_ch)] = aln
            return MappingResult(self.target, mdl, self.chem_groups, chem_mapping,
                                 unmapped_mdl_chains, one_to_one, alns)

        if strategy == "heuristic":
            if _NMappingsWithin(self.chem_groups, chem_mapping,
                                heuristic_n_max_naive):
                strategy = "naive"
            else:
                strategy = "greedy_full"

        mapping = None
        opt_lddt = None

        if strategy == "naive":
            mapping, opt_lddt = _lDDTNaive(self.target, mdl, inclusion_radius,
                                           thresholds, self.chem_groups,
                                           chem_mapping, ref_mdl_alns,
                                           self.n_max_naive)
        else:
            # its one of the greedy strategies - setup greedy searcher
            the_greed = _lDDTGreedySearcher(self.target, mdl, self.chem_groups,
                                            chem_mapping, ref_mdl_alns,
                                            inclusion_radius=inclusion_radius,
                                            thresholds=thresholds,
                                            steep_opt_rate=steep_opt_rate)
            if strategy == "greedy_fast":
                mapping = _lDDTGreedyFast(the_greed)
            elif strategy == "greedy_full":
                mapping = _lDDTGreedyFull(the_greed)
            elif strategy == "greedy_block":
                mapping = _lDDTGreedyBlock(the_greed, block_seed_size,
                                           block_blocks_per_chem_group)
            # cached => lDDT computation is fast here
            opt_lddt = the_greed.Score(mapping)

        alns = dict()
        for ref_group, mdl_group in zip(self.chem_groups, mapping):
            for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                if ref_ch is not None and mdl_ch is not None:
                    aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                    alns[(ref_ch, mdl_ch)] = aln

        return MappingResult(self.target, mdl, self.chem_groups, chem_mapping,
                             unmapped_mdl_chains, mapping, alns,
                             opt_score = opt_lddt)


    def GetQSScoreMapping(self, model, contact_d = 12.0, strategy = "heuristic",
                          block_seed_size = 5, block_blocks_per_chem_group = 5,
                          heuristic_n_max_naive = 40320, steep_opt_rate = None,
                          chem_mapping_result = None,
                          greedy_prune_contact_map = True):
        """ Identify chain mapping based on QSScore

        Scoring is based on CA/C3' positions which are present in all chains of
        a :attr:`chem_groups` as well as the *model* chains which are mapped to
        that respective chem group.

        The following strategies are available:

        * **naive**: Naively iterate all possible mappings and return best based
                     on QS score.

        * **greedy_fast**: perform all vs. all single chain lDDTs within the
          respective ref/mdl chem groups. The mapping with highest number of
          conserved contacts is selected as seed for greedy extension.
          Extension is based on QS-score.

        * **greedy_full**: try multiple seeds for greedy extension, i.e. try
          all ref/mdl chain combinations within the respective chem groups and
          retain the mapping leading to the best QS-score. 

        * **greedy_block**: try multiple seeds for greedy extension, i.e. try
          all ref/mdl chain combinations within the respective chem groups and
          extend them to *block_seed_size*. *block_blocks_per_chem_group*
          for each chem group are selected for exhaustive extension.

        * **heuristic**: Uses *naive* strategy if number of possible mappings
          is within *heuristic_n_max_naive*. The default of 40320 corresponds
          to an octamer (8!=40320). A structure with stoichiometry A6B2 would be
          6!*2!=1440 etc. If the number of possible mappings is larger,
          *greedy_full* is used.

        Sets :attr:`MappingResult.opt_score` in case of no trivial one-to-one
        mapping.

        :param model: Model to map
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param contact_d: Max distance between two residues to be considered as 
                          contact in qs scoring
        :type contact_d: :class:`float` 
        :param strategy: Strategy for sampling, must be in ["naive",
                         "greedy_fast", "greedy_full", "greedy_block"]
        :type strategy: :class:`str`
        :param chem_mapping_result: Pro param. The result of
                                    :func:`~GetChemMapping` where you provided
                                    *model*. If set, *model* parameter is not
                                    used.
        :type chem_mapping_result: :class:`tuple`
        :param greedy_prune_contact_map: Relevant for all strategies that use
                                         greedy extensions. If True, only chains
                                         with at least 3 contacts (8A CB
                                         distance) towards already mapped chains
                                         in trg/mdl are considered for
                                         extension. All chains that give a
                                         potential non-zero QS-score increase
                                         are used otherwise (at least one
                                         contact within 12A). The consequence
                                         is reduced runtime and usually no
                                         real reduction in accuracy.
        :returns: A :class:`MappingResult`
        """

        strategies = ["naive", "greedy_fast", "greedy_full", "greedy_block", "heuristic"]
        if strategy not in strategies:
            raise RuntimeError(f"strategy must be {strategies}")

        if chem_mapping_result is None:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            self.GetChemMapping(model)
        else:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            chem_mapping_result
        ref_mdl_alns =  _GetRefMdlAlns(self.chem_groups,
                                       self.chem_group_alignments,
                                       chem_mapping,
                                       chem_group_alns)
        # check for the simplest case
        one_to_one = _CheckOneToOneMapping(self.chem_groups, chem_mapping)
        if one_to_one is not None:
            alns = dict()
            for ref_group, mdl_group in zip(self.chem_groups, one_to_one):
                for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                    if ref_ch is not None and mdl_ch is not None:
                        aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                        alns[(ref_ch, mdl_ch)] = aln
            return MappingResult(self.target, mdl, self.chem_groups, chem_mapping,
                                 unmapped_mdl_chains, one_to_one, alns)

        if strategy == "heuristic":
            if _NMappingsWithin(self.chem_groups, chem_mapping,
                                heuristic_n_max_naive):
                strategy = "naive"
            else:
                strategy = "greedy_full"

        mapping = None
        opt_qsscore = None

        if strategy == "naive":
            mapping, opt_qsscore = _QSScoreNaive(self.target, mdl,
                                                 self.chem_groups,
                                                 chem_mapping, ref_mdl_alns,
                                                 contact_d, self.n_max_naive)
        else:
            # its one of the greedy strategies - setup greedy searcher
            the_greed = _QSScoreGreedySearcher(self.target, mdl,
                                               self.chem_groups,
                                               chem_mapping, ref_mdl_alns,
                                               contact_d = contact_d,
                                               steep_opt_rate=steep_opt_rate,
                                               greedy_prune_contact_map = greedy_prune_contact_map)
            if strategy == "greedy_fast":
                mapping = _QSScoreGreedyFast(the_greed)
            elif strategy == "greedy_full":
                mapping = _QSScoreGreedyFull(the_greed)
            elif strategy == "greedy_block":
                mapping = _QSScoreGreedyBlock(the_greed, block_seed_size,
                                              block_blocks_per_chem_group)
            # cached => QSScore computation is fast here
            opt_qsscore = the_greed.Score(mapping, check=False).QS_global
              
        alns = dict()
        for ref_group, mdl_group in zip(self.chem_groups, mapping):
            for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                if ref_ch is not None and mdl_ch is not None:
                    aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                    alns[(ref_ch, mdl_ch)] = aln

        return MappingResult(self.target, mdl, self.chem_groups, chem_mapping,
                             unmapped_mdl_chains, mapping, alns,
                             opt_score = opt_qsscore)

    def GetRMSDMapping(self, model, strategy = "heuristic", subsampling=50,
                       chem_mapping_result = None, heuristic_n_max_naive = 120):
        """Identify chain mapping based on minimal RMSD superposition

        Superposition and scoring is based on CA/C3' positions which are present
        in all chains of a :attr:`chem_groups` as well as the *model*
        chains which are mapped to that respective chem group.

        The following strategies are available:

        * **naive**: Naively iterate all possible mappings and return the one
          with lowes RMSD.

        * **greedy_single**: perform all vs. all single chain superpositions
          within the respective ref/mdl chem groups to use as starting points.
          For each starting point, iteratively add the model/target chain pair
          with lowest RMSD until a full mapping is achieved. The mapping with
          lowest RMSD is returned.

        * **greedy_iterative**: Same as greedy_single_rmsd exept that the
          transformation gets updated with each added chain pair.

        * **heuristic**: Uses *naive* strategy if number of possible mappings
          is within *heuristic_n_max_naive*. The default of 120 corresponds
          to a homo-pentamer (5!=120). If the number of possible mappings is
          larger, *greedy_iterative* is used.

        :param model: Model to map
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param strategy: Strategy for sampling. Must be in ["naive",
                         "greedy_single", "greedy_iterative"]
        :type strategy: :class:`str`
        :param subsampling: If given, only an equally distributed subset
                            of CA/C3' positions in each chain are used for
                            superposition/scoring.
        :type subsampling: :class:`int`
        :param chem_mapping_result: Pro param. The result of
                                    :func:`~GetChemMapping` where you provided
                                    *model*. If set, *model* parameter is not
                                    used.
        :type chem_mapping_result: :class:`tuple`
        :returns: A :class:`MappingResult`
        """

        strategies = ["naive", "greedy_single", "greedy_iterative", "heuristic"]

        if strategy not in strategies:
            raise RuntimeError(f"strategy must be {strategies}")

        if chem_mapping_result is None:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            self.GetChemMapping(model)
        else:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            chem_mapping_result
        ref_mdl_alns =  _GetRefMdlAlns(self.chem_groups,
                                       self.chem_group_alignments,
                                       chem_mapping,
                                       chem_group_alns)

        # check for the simplest case
        one_to_one = _CheckOneToOneMapping(self.chem_groups, chem_mapping)
        if one_to_one is not None:
            alns = dict()
            for ref_group, mdl_group in zip(self.chem_groups, one_to_one):
                for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                    if ref_ch is not None and mdl_ch is not None:
                        aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                        alns[(ref_ch, mdl_ch)] = aln
            return MappingResult(self.target, mdl, self.chem_groups, chem_mapping,
                                 unmapped_mdl_chains, one_to_one, alns)

        trg_group_pos, mdl_group_pos = _GetRefPos(self.target, mdl,
                                                  self.chem_group_alignments,
                                                  chem_group_alns,
                                                  max_pos = subsampling)

        if strategy == "heuristic":
            if _NMappingsWithin(self.chem_groups, chem_mapping,
                                heuristic_n_max_naive):
                strategy = "naive"
            else:
                strategy = "greedy_iterative"

        mapping = None

        if strategy.startswith("greedy"):
          # get transforms of any mdl chain onto any trg chain in same chem
          # group that fulfills gdtts threshold
          initial_transforms = list()
          initial_mappings = list()
          for trg_pos, trg_chains, mdl_pos, mdl_chains in zip(trg_group_pos,
                                                              self.chem_groups,
                                                              mdl_group_pos,
                                                              chem_mapping):
              for t_pos, t in zip(trg_pos, trg_chains):
                  for m_pos, m in zip(mdl_pos, mdl_chains):
                      if len(t_pos) >= 3 and len(m_pos) >= 3:
                          transform = _GetSuperposition(m_pos, t_pos,
                                                        False).transformation
                          initial_transforms.append(transform)
                          initial_mappings.append((t,m))

          if strategy == "greedy_single":
              mapping = _SingleRigidRMSD(initial_transforms,
                                         initial_mappings,
                                         self.chem_groups,
                                         chem_mapping,
                                         trg_group_pos,
                                         mdl_group_pos)


          elif strategy == "greedy_iterative":
              mapping = _IterativeRigidRMSD(initial_transforms,
                                            initial_mappings,
                                            self.chem_groups, chem_mapping,
                                            trg_group_pos, mdl_group_pos)
        elif strategy == "naive":
            mapping = _NaiveRMSD(self.chem_groups, chem_mapping,
                                 trg_group_pos, mdl_group_pos,
                                 self.n_max_naive)

        # translate mapping format and return
        final_mapping = list()
        for ref_chains in self.chem_groups:
            mapped_mdl_chains = list()
            for ref_ch in ref_chains:
                if ref_ch in mapping:
                    mapped_mdl_chains.append(mapping[ref_ch])
                else:
                    mapped_mdl_chains.append(None)
            final_mapping.append(mapped_mdl_chains)

        alns = dict()
        for ref_group, mdl_group in zip(self.chem_groups, final_mapping):
            for ref_ch, mdl_ch in zip(ref_group, mdl_group):
                if ref_ch is not None and mdl_ch is not None:
                    aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                    alns[(ref_ch, mdl_ch)] = aln

        return MappingResult(self.target, mdl, self.chem_groups, chem_mapping,
                             unmapped_mdl_chains, final_mapping, alns)

    def GetMapping(self, model, n_max_naive = 40320):
        """ Convenience function to get mapping with currently preferred method

        If number of possible chain mappings is <= *n_max_naive*, a naive
        QS-score mapping is performed and optimal QS-score is guaranteed.
        For anything else, a QS-score mapping with the greedy_full strategy is
        performed (greedy_prune_contact_map = True). The default for
        *n_max_naive* of 40320 corresponds to an octamer (8!=40320). A
        structure with stoichiometry A6B2 would be 6!*2!=1440 etc.

        If :attr:`~target` has nucleotide sequences, the QS-score target
        function is replaced with a backbone only lDDT score that has
        an inclusion radius of 30A. 
        """
        if len(self.polynuc_seqs) > 0:
            return self.GetlDDTMapping(model, strategy = "heuristic",
                                       inclusion_radius = 30.0,
                                       heuristic_n_max_naive = n_max_naive)
        else:
            return self.GetQSScoreMapping(model, strategy="heuristic",
                                          greedy_prune_contact_map=True,
                                          heuristic_n_max_naive = n_max_naive)

    def GetRepr(self, substructure, model, topn=1, inclusion_radius=15.0,
                thresholds=[0.5, 1.0, 2.0, 4.0], bb_only=False,
                only_interchain=False, chem_mapping_result = None,
                global_mapping = None):
        """ Identify *topn* representations of *substructure* in *model*

        *substructure* defines a subset of :attr:`~target` for which one
        wants the *topn* representations in *model*. Representations are scored
        and sorted by lDDT.

        :param substructure: A :class:`ost.mol.EntityView` which is a subset of
                             :attr:`~target`. Should be selected with the
                             OpenStructure query language. Example: if you're
                             interested in residues with number 42,43 and 85 in
                             chain A:
                             ``substructure=mapper.target.Select("cname=A and rnum=42,43,85")``
                             A :class:`RuntimeError` is raised if *substructure*
                             does not refer to the same underlying
                             :class:`ost.mol.EntityHandle` as :attr:`~target`.
        :type substructure: :class:`ost.mol.EntityView`
        :param model: Structure in which one wants to find representations for
                      *substructure*
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param topn: Max number of representations that are returned
        :type topn: :class:`int`
        :param inclusion_radius: Inclusion radius for lDDT
        :type inclusion_radius: :class:`float`
        :param thresholds: Thresholds for lDDT
        :type thresholds: :class:`list` of :class:`float`
        :param bb_only: Only consider backbone atoms in lDDT computation
        :type bb_only: :class:`bool`
        :param only_interchain: Only score interchain contacts in lDDT. Useful
                                if you want to identify interface patches.
        :type only_interchain: :class:`bool`
        :param chem_mapping_result: Pro param. The result of
                                    :func:`~GetChemMapping` where you provided
                                    *model*. If set, *model* parameter is not
                                    used.
        :type chem_mapping_result: :class:`tuple`
        :param global_mapping: Pro param. Specify a global mapping result. This
                               fully defines the desired representation in the
                               model but extracts it and enriches it with all
                               the nice attributes of :class:`ReprResult`.
                               The target attribute in *global_mapping* must be
                               of the same entity as self.target and the model
                               attribute of *global_mapping* must be of the same
                               entity as *model*.
        :type global_mapping: :class:`MappingResult`
        :returns: :class:`list` of :class:`ReprResult`
        """

        if topn < 1:
            raise RuntimeError("topn must be >= 1")

        if global_mapping is not None:
            # ensure that this mapping is derived from the same structures
            if global_mapping.target.handle.GetHashCode() != \
               self.target.handle.GetHashCode():
               raise RuntimeError("global_mapping.target must be the same "
                                  "entity as self.target")
            if global_mapping.model.handle.GetHashCode() != \
               model.handle.GetHashCode():
               raise RuntimeError("global_mapping.model must be the same "
                                  "entity as model param")

        # check whether substructure really is a subset of self.target
        for r in substructure.residues:
            ch_name = r.GetChain().GetName()
            rnum = r.GetNumber()
            target_r = self.target.FindResidue(ch_name, rnum)
            if not target_r.IsValid():
                raise RuntimeError(f"substructure has residue "
                                   f"{r.GetQualifiedName()} which is not in "
                                   f"self.target")
            if target_r.handle.GetHashCode() != r.handle.GetHashCode():
                raise RuntimeError(f"substructure has residue "
                                   f"{r.GetQualifiedName()} which has an "
                                   f"equivalent in self.target but it does "
                                   f"not refer to the same underlying "
                                   f"EntityHandle")
            for a in r.atoms:
                target_a = target_r.FindAtom(a.GetName())
                if not target_a.IsValid():
                    raise RuntimeError(f"substructure has atom "
                                       f"{a.GetQualifiedName()} which is not "
                                       f"in self.target")
                if a.handle.GetHashCode() != target_a.handle.GetHashCode():
                    raise RuntimeError(f"substructure has atom "
                                       f"{a.GetQualifiedName()} which has an "
                                       f"equivalent in self.target but it does "
                                       f"not refer to the same underlying "
                                       f"EntityHandle")

            # check whether it contains either CA or C3'
            ca = r.FindAtom("CA")
            c3 = r.FindAtom("C3'") # FindAtom with prime in string is tested
                                   # and works
            if not ca.IsValid() and not c3.IsValid():
                raise RuntimeError("All residues in substructure must contain "
                                   "a backbone atom named CA or C3\'")

        # perform mapping and alignments on full structures
        if chem_mapping_result is None:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            self.GetChemMapping(model)
        else:
            chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
            chem_mapping_result
        ref_mdl_alns =  _GetRefMdlAlns(self.chem_groups,
                                       self.chem_group_alignments,
                                       chem_mapping,
                                       chem_group_alns)

        # Get residue indices relative to full target chain 
        substructure_res_indices = dict()
        for ch in substructure.chains:
            full_ch = self.target.FindChain(ch.GetName())
            idx = [full_ch.GetResidueIndex(r.GetNumber()) for r in ch.residues]
            substructure_res_indices[ch.GetName()] = idx

        # strip down variables to make them specific to substructure
        # keep only chem_groups which are present in substructure
        substructure_chem_groups = list()
        substructure_chem_mapping = list()
        
        chnames = set([ch.GetName() for ch in substructure.chains])
        for chem_group, mapping in zip(self.chem_groups, chem_mapping):
            substructure_chem_group = [ch for ch in chem_group if ch in chnames]
            if len(substructure_chem_group) > 0:
                substructure_chem_groups.append(substructure_chem_group)
                substructure_chem_mapping.append(mapping)

        # early stopping if no mdl chain can be mapped to substructure
        n_mapped_mdl_chains = sum([len(m) for m in substructure_chem_mapping])
        if n_mapped_mdl_chains == 0:
            return list()

        # strip the reference sequence in alignments to only contain
        # sequence from substructure
        substructure_ref_mdl_alns = dict()
        mdl_views = dict()
        for ch in mdl.chains:
            mdl_views[ch.GetName()] = _CSel(mdl, [ch.GetName()])
        for chem_group, mapping in zip(substructure_chem_groups,
                                       substructure_chem_mapping):
            for ref_ch in chem_group:
                for mdl_ch in mapping:
                    full_aln = ref_mdl_alns[(ref_ch, mdl_ch)]
                    ref_seq = full_aln.GetSequence(0)
                    # the ref sequence is tricky... we start with a gap only
                    # sequence and only add olcs as defined by the residue
                    # indices that we extracted before...
                    tmp = ['-'] * len(full_aln)
                    for idx in substructure_res_indices[ref_ch]:
                        idx_in_seq = ref_seq.GetPos(idx)
                        tmp[idx_in_seq] = ref_seq[idx_in_seq]
                    ref_seq = seq.CreateSequence(ref_ch, ''.join(tmp))
                    ref_seq.AttachView(_CSel(substructure, [ref_ch]))
                    mdl_seq = full_aln.GetSequence(1)
                    mdl_seq = seq.CreateSequence(mdl_seq.GetName(),
                                                 mdl_seq.GetString())
                    mdl_seq.AttachView(mdl_views[mdl_ch])
                    aln = seq.CreateAlignment()
                    aln.AddSequence(ref_seq)
                    aln.AddSequence(mdl_seq)
                    substructure_ref_mdl_alns[(ref_ch, mdl_ch)] = aln

        lddt_scorer = lddt.lDDTScorer(substructure,
                                      inclusion_radius = inclusion_radius,
                                      bb_only = bb_only)
        scored_mappings = list()

        if global_mapping:
            # construct mapping of substructure from global mapping
            flat_mapping = global_mapping.GetFlatMapping()
            mapping = list()
            for chem_group, chem_mapping in zip(substructure_chem_groups,
                                                substructure_chem_mapping):
                chem_group_mapping = list()
                for ch in chem_group:
                    if ch in flat_mapping:
                        mdl_ch = flat_mapping[ch]
                        if mdl_ch in chem_mapping:
                            chem_group_mapping.append(mdl_ch)
                        else:
                            chem_group_mapping.append(None)
                    else:
                        chem_group_mapping.append(None)
                mapping.append(chem_group_mapping)
            mappings = [mapping]
        else:
            mappings = list(_ChainMappings(substructure_chem_groups,
                                           substructure_chem_mapping,
                                           self.n_max_naive))

        # This step can be slow so give some hints in logs
        msg = "Computing initial ranking of %d chain mappings" % len(mappings)
        (ost.LogWarning if len(mappings) > 10000 else ost.LogInfo)(msg)

        for mapping in mappings:
            # chain_mapping and alns as input for lDDT computation
            lddt_chain_mapping = dict()
            lddt_alns = dict()
            n_res_aln = 0
            for ref_chem_group, mdl_chem_group in zip(substructure_chem_groups,
                                                      mapping):
                for ref_ch, mdl_ch in zip(ref_chem_group, mdl_chem_group):
                    # some mdl chains can be None
                    if mdl_ch is not None:
                        lddt_chain_mapping[mdl_ch] = ref_ch
                        aln = substructure_ref_mdl_alns[(ref_ch, mdl_ch)]
                        lddt_alns[mdl_ch] = aln
                        tmp = [int(c[0] != '-' and c[1] != '-') for c in aln]
                        n_res_aln += sum(tmp)
            # don't compute lDDT if no single residue in mdl and ref is aligned
            if n_res_aln == 0:
                continue

            lDDT, _ = lddt_scorer.lDDT(mdl, thresholds=thresholds,
                                       chain_mapping=lddt_chain_mapping,
                                       residue_mapping = lddt_alns,
                                       check_resnames = False,
                                       no_intrachain = only_interchain)

            if lDDT is None:
                ost.LogVerbose("No valid contacts in the reference")
                lDDT = 0.0 # that means, that we have not a single valid contact
                           # in lDDT. For the code below to work, we just set it
                           # to a terrible score => 0.0

            if len(scored_mappings) == 0:
                scored_mappings.append((lDDT, mapping))
            elif len(scored_mappings) < topn:
                scored_mappings.append((lDDT, mapping))
                scored_mappings.sort(reverse=True, key=lambda x: x[0])
            elif lDDT > scored_mappings[-1][0]:
                scored_mappings.append((lDDT, mapping))
                scored_mappings.sort(reverse=True, key=lambda x: x[0])
                scored_mappings = scored_mappings[:topn]

        # finalize and return
        results = list()
        for scored_mapping in scored_mappings:
            ref_view = substructure.handle.CreateEmptyView()
            mdl_view = mdl.handle.CreateEmptyView()
            for ref_ch_group, mdl_ch_group in zip(substructure_chem_groups,
                                                  scored_mapping[1]):
                for ref_ch, mdl_ch in zip(ref_ch_group, mdl_ch_group):
                    if ref_ch is not None and mdl_ch is not None:
                        aln = substructure_ref_mdl_alns[(ref_ch, mdl_ch)]
                        for col in aln:
                            if col[0] != '-' and col[1] != '-':
                                ref_view.AddResidue(col.GetResidue(0),
                                                    mol.ViewAddFlag.INCLUDE_ALL)
                                mdl_view.AddResidue(col.GetResidue(1),
                                                    mol.ViewAddFlag.INCLUDE_ALL)
            results.append(ReprResult(scored_mapping[0], substructure,
                                      ref_view, mdl_view))
        return results

    def GetNMappings(self, model):
        """ Returns number of possible mappings

        :param model: Model with chains that are mapped onto
                      :attr:`chem_groups`
        :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        """
        chem_mapping, chem_group_alns, unmapped_mdl_chains, mdl = \
        self.GetChemMapping(model)
        return _NMappings(self.chem_groups, chem_mapping)

    def ProcessStructure(self, ent):
        """ Entity processing for chain mapping

        * Selects view containing peptide and nucleotide residues which have 
          required backbone atoms present - for peptide residues thats
          N, CA, C and CB (no CB for GLY), for nucleotide residues thats
          O5', C5', C4', C3' and O3'.
        * filters view by chain lengths, see *min_pep_length* and
          *min_nuc_length* in constructor
        * Extracts atom sequences for each chain in that view
        * If residue number alignments are used, strictly increasing residue
          numbers without insertion codes are ensured in each chain

        :param ent: Entity to process
        :type ent: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :returns: Tuple with 3 elements: 1) :class:`ost.mol.EntityView`
                  containing peptide and nucleotide residues 2)
                  :class:`ost.seq.SequenceList` containing ATOMSEQ sequences
                  for each polypeptide chain in returned view
                  3) same for polynucleotide chains
        """
        view = ent.CreateEmptyView()
        exp_pep_atoms = ["N", "CA", "C", "CB"]
        exp_nuc_atoms = ["\"O5'\"", "\"C5'\"", "\"C4'\"", "\"C3'\"", "\"O3'\""]
        pep_query = "peptide=true and aname=" + ','.join(exp_pep_atoms)
        nuc_query = "nucleotide=true and aname=" + ','.join(exp_nuc_atoms)

        pep_sel = ent.Select(pep_query)
        for r in pep_sel.residues:
            if len(r.atoms) == 4:
                view.AddResidue(r.handle, mol.INCLUDE_ALL)
            elif r.name == "GLY" and len(r.atoms) == 3:
                atom_names = [a.GetName() for a in r.atoms]
                if sorted(atom_names) == ["C", "CA", "N"]:
                    view.AddResidue(r.handle, mol.INCLUDE_ALL)

        nuc_sel = ent.Select(nuc_query)
        for r in nuc_sel.residues:
            if len(r.atoms) == 5:
                view.AddResidue(r.handle, mol.INCLUDE_ALL)

        polypep_seqs = seq.CreateSequenceList()
        polynuc_seqs = seq.CreateSequenceList()

        if len(view.residues) == 0:
            # no residues survived => return
            return (view, polypep_seqs, polynuc_seqs)

        for ch in view.chains:
            n_res = len(ch.residues)
            n_pep = sum([r.IsPeptideLinking() for r in ch.residues])
            n_nuc = sum([r.IsNucleotideLinking() for r in ch.residues])

            # guarantee that we have either pep or nuc (no mix of the two)
            if n_pep > 0 and n_nuc > 0:
                raise RuntimeError(f"Must not mix peptide and nucleotide linking "
                                   f"residues in same chain ({ch.GetName()})")

            if (n_pep + n_nuc) != n_res:
                raise RuntimeError("All residues must either be peptide_linking "
                                   "or nucleotide_linking")

            # filter out short chains
            if n_pep > 0 and n_pep < self.min_pep_length:
                continue

            if n_nuc > 0 and n_nuc < self.min_nuc_length:
                continue

            # the superfast residue number based alignment adds some 
            # restrictions on the numbers themselves:
            # 1) no insertion codes 2) strictly increasing
            if self.resnum_alignments:
                # check if no insertion codes are present in residue numbers
                ins_codes = [r.GetNumber().GetInsCode() for r in ch.residues]
                if len(set(ins_codes)) != 1 or ins_codes[0] != '\0':
                    raise RuntimeError("Residue numbers in input structures must not "
                                       "contain insertion codes")

                # check if residue numbers are strictly increasing
                nums = [r.GetNumber().GetNum() for r in ch.residues]
                if not all(i < j for i, j in zip(nums, nums[1:])):
                    raise RuntimeError("Residue numbers in input structures must be "
                                       "strictly increasing for each chain")

            s = ''.join([r.one_letter_code for r in ch.residues])
            s = seq.CreateSequence(ch.GetName(), s)
            if n_pep == n_res:
                polypep_seqs.AddSequence(s)
            elif n_nuc == n_res:
                polynuc_seqs.AddSequence(s)
            else:
                raise RuntimeError("This shouldnt happen")

        if len(polypep_seqs) == 0 and len(polynuc_seqs) == 0:
            raise RuntimeError(f"No chain fulfilled minimum length requirement "
                               f"to be considered in chain mapping "
                               f"({self.min_pep_length} for peptide chains, "
                               f"{self.min_nuc_length} for nucleotide chains) "
                               f"- mapping failed")

        # select for chains for which we actually extracted the sequence
        chain_names = [s.name for s in polypep_seqs]
        chain_names += [s.name for s in polynuc_seqs]
        view = _CSel(view, chain_names)

        return (view, polypep_seqs, polynuc_seqs)

    def NWAlign(self, s1, s2, stype):
        """ Access to internal sequence alignment functionality

        Performs Needleman-Wunsch alignment with parameterization
        setup at ChainMapper construction

        :param s1: First sequence to align
        :type s1: :class:`ost.seq.SequenceHandle`
        :param s2: Second sequence to align
        :type s2: :class:`ost.seq.SequenceHandle`
        :param stype: Type of sequences to align, must be in
                      [:class:`ost.mol.ChemType.AMINOACIDS`,
                      :class:`ost.mol.ChemType.NUCLEOTIDES`]
        :returns: Pairwise alignment of s1 and s2
        """
        if stype == mol.ChemType.AMINOACIDS:
            return seq.alg.SemiGlobalAlign(s1, s2, self.pep_subst_mat,
                                           gap_open=self.pep_gap_open,
                                           gap_ext=self.pep_gap_ext)[0]
        elif stype == mol.ChemType.NUCLEOTIDES:
            return seq.alg.SemiGlobalAlign(s1, s2, self.nuc_subst_mat,
                                           gap_open=self.nuc_gap_open,
                                           gap_ext=self.nuc_gap_ext)[0]
        else:
            raise RuntimeError("Invalid ChemType")
        return aln

    def ResNumAlign(self, s1, s2, s1_ent, s2_ent):
        """ Access to internal sequence alignment functionality

        Performs residue number based alignment. Residue numbers are extracted
        from *s1_ent*/*s2_ent*.

        :param s1: First sequence to align
        :type s1: :class:`ost.seq.SequenceHandle`
        :param s2: Second sequence to align
        :type s2: :class:`ost.seq.SequenceHandle`
        :param s1_ent: Structure as source for residue numbers in *s1*. Must
                       contain a chain named after sequence name in *s1*. This
                       chain must have the exact same number of residues as
                       characters in s1.
        :type s1_ent: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
        :param s2_ent: Same for *s2*.
        :type s2_ent: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`

        :returns: Pairwise alignment of s1 and s2
        """
        ch1 = s1_ent.FindChain(s1.name)
        ch2 = s2_ent.FindChain(s2.name)

        if not ch1.IsValid():
            raise RuntimeError("s1_ent lacks requested chain in ResNumAlign")

        if not ch2.IsValid():
            raise RuntimeError("s2_ent lacks requested chain in ResNumAlign")

        if len(s1) != ch1.GetResidueCount():
            raise RuntimeError("Sequence/structure mismatch in ResNumAlign")

        if len(s2) != ch2.GetResidueCount():
            raise RuntimeError("Sequence/structure mismatch in ResNumAlign")

        rnums1 = [r.GetNumber().GetNum() for r in ch1.residues]
        rnums2 = [r.GetNumber().GetNum() for r in ch2.residues]

        min_num = min(rnums1[0], rnums2[0])
        max_num = max(rnums1[-1], rnums2[-1])
        aln_length = max_num - min_num + 1

        aln_s1 = ['-'] * aln_length
        for olc, rnum in zip(s1, rnums1):
            aln_s1[rnum-min_num] = olc

        aln_s2 = ['-'] * aln_length
        for olc, rnum in zip(s2, rnums2):
            aln_s2[rnum-min_num] = olc

        aln = seq.CreateAlignment()
        aln.AddSequence(seq.CreateSequence(s1.GetName(), ''.join(aln_s1)))
        aln.AddSequence(seq.CreateSequence(s2.GetName(), ''.join(aln_s2)))
        return aln


    def _ComputeChemGroups(self):
        """ Sets properties :attr:`~chem_groups`,
        :attr:`~chem_group_alignments`, :attr:`~chem_group_ref_seqs`,
        :attr:`~chem_group_types`
        """

        self._chem_group_alignments, self._chem_group_types =\
        self._ChemGroupAlignmentsFromATOMSEQ()


        self._chem_group_ref_seqs = seq.CreateSequenceList()
        for a in self.chem_group_alignments:
            s = seq.CreateSequence(a.GetSequence(0).GetName(),
                                   a.GetSequence(0).GetGaplessString())
            self._chem_group_ref_seqs.AddSequence(s)

        self._chem_groups = list()
        for a in self.chem_group_alignments:
            group = list()
            for s_idx in range(1, a.GetCount()):
                s = a.GetSequence(s_idx)
                group.append(s.GetName())
            self._chem_groups.append(group)

    def _ChemGroupAlignmentsFromATOMSEQ(self):
        """ Groups target sequences based on ATOMSEQ

        returns tuple that can be set as self._chem_group_alignments and
        self._chem_group_types
        """
        pep_groups = self._GroupSequences(self.polypep_seqs, self.pep_seqid_thr,
                                          self.min_pep_length,
                                          mol.ChemType.AMINOACIDS)
        nuc_groups = self._GroupSequences(self.polynuc_seqs, self.nuc_seqid_thr,
                                          self.min_nuc_length,
                                          mol.ChemType.NUCLEOTIDES)

        # pep_groups and nuc_groups give us alignments based on ATOMSEQ.
        # For example: If you have polymer chains A,B and C in the same
        # group and A is the longest one, you get an alignment that looks
        # like:
        # A: ASDFE
        # B: ASDF-
        # C: -SDFE
        #
        # however, the first sequence in chem group alignments must not be
        # bound to any ATOMSEQ and represent the reference sequence. In the
        # case of this function, this is simply a copy of sequence A:
        # REF: ASDFE
        # A:   ASDFE
        # B:   ASDF-
        # C:   -SDFE

        # do pep_groups
        tmp = list()
        for a in pep_groups:
            new_a = seq.CreateAlignment()
            new_a.AddSequence(a.GetSequence(0))
            for s_idx in range(a.GetCount()):
                new_a.AddSequence(a.GetSequence(s_idx))
            tmp.append(new_a)
        pep_groups = tmp

        # do nuc groups        
        tmp = list()
        for a in nuc_groups:
            new_a = seq.CreateAlignment()
            new_a.AddSequence(a.GetSequence(0))
            for s_idx in range(a.GetCount()):
                new_a.AddSequence(a.GetSequence(s_idx))
            tmp.append(new_a)
        nuc_groups = tmp

        group_types = [mol.ChemType.AMINOACIDS] * len(pep_groups)
        group_types += [mol.ChemType.NUCLEOTIDES] * len(nuc_groups)
        groups = pep_groups
        groups.extend(nuc_groups)

        return (groups, group_types)

    def _GroupSequences(self, seqs, seqid_thr, min_length, chem_type):
        """Get list of alignments representing groups of equivalent sequences
    
        :param seqid_thr: Threshold used to decide when two chains are identical.
        :type seqid_thr:  :class:`float`
        :param gap_thr: Additional threshold to avoid gappy alignments with high
                        seqid. Number of aligned columns must be at least this
                        number.
        :type gap_thr: :class:`int`
        :param aligner: Helper class to generate pairwise alignments
        :type aligner: :class:`_Aligner`
        :param chem_type: ChemType of seqs, must be in
                          [:class:`ost.mol.ChemType.AMINOACIDS`,
                          :class:`ost.mol.ChemType.NUCLEOTIDES`]
        :type chem_type: :class:`ost.mol.ChemType` 
        :returns: A list of alignments, one alignment for each group
                  with longest sequence (reference) as first sequence.
        :rtype: :class:`ost.seq.AlignmentList`
        """
        groups = list()
        for s_idx in range(len(seqs)):
            matching_group = None
            for g_idx in range(len(groups)):
                for g_s_idx in range(len(groups[g_idx])):
                    if self.resnum_alignments:
                        aln = self.ResNumAlign(seqs[s_idx],
                                               seqs[groups[g_idx][g_s_idx]],
                                               self.target, self.target)
                    else:
                        aln = self.NWAlign(seqs[s_idx],
                                           seqs[groups[g_idx][g_s_idx]],
                                           chem_type)
                    sid, n_aligned = _GetAlnPropsOne(aln)
                    if sid >= seqid_thr and n_aligned >= min_length:
                        matching_group = g_idx
                        break
                if matching_group is not None:
                    break
    
            if matching_group is None:
                groups.append([s_idx])
            else:
                groups[matching_group].append(s_idx)
    
        # sort based on sequence length
        sorted_groups = list()
        for g in groups:
            if len(g) > 1:
                tmp = sorted([[len(seqs[i]), i] for i in g], reverse=True)
                sorted_groups.append([x[1] for x in tmp])
            else:
                sorted_groups.append(g)
    
        # translate from indices back to sequences and directly generate alignments
        # of the groups with the longest (first) sequence as reference
        aln_list = seq.AlignmentList()
        for g in sorted_groups:
            if len(g) == 1:
                # aln with one single sequence
                aln_list.append(seq.CreateAlignment(seqs[g[0]]))
            else:
                # obtain pairwise aln of first sequence (reference) to all others
                alns = seq.AlignmentList()
                i = g[0]
                for j in g[1:]:
                    if self.resnum_alignments:
                        aln = self.ResNumAlign(seqs[i], seqs[j],
                                               self.target, self.target)
                    else:
                        aln = self.NWAlign(seqs[i], seqs[j], chem_type)
                    alns.append(aln)
                # and merge
                aln_list.append(seq.alg.MergePairwiseAlignments(alns, seqs[i]))
    
        return aln_list

    def _MapSequence(self, ref_seqs, ref_types, s, s_type, s_ent,
                     seq_id_thr=0.0, min_aln_length=0):
        """Tries top map *s* onto any of the sequences in *ref_seqs*
    
        Computes alignments of *s* to each of the reference sequences of equal type
        and sorts them by seqid*fraction_covered (seqid: sequence identity of
        aligned columns in alignment, fraction_covered: Fraction of non-gap
        characters in reference sequence that are covered by non-gap characters in
        *s*). Best scoring mapping is returned. Optionally, *seq_id*/
        *min_aln_length* thresholds can be enabled to avoid non-sensical mappings.
        However, *min_aln_length* only has an effect if *seq_id_thr* > 0!!!
    
        :param ref_seqs: Reference sequences 
        :type ref_seqs: :class:`ost.seq.SequenceList`
        :param ref_types: Types of reference sequences, e.g.
                          ost.mol.ChemType.AminoAcids
        :type ref_types: :class:`list` of :class:`ost.mol.ChemType`
        :param s: Sequence to map
        :type s: :class:`ost.seq.SequenceHandle`
        :param s_type: Type of *s*, only try mapping to sequences in *ref_seqs*
                       with equal type as defined in *ref_types*
        :param s_ent: Entity which represents *s*. Only relevant in case of
                      residue number alignments.
        :type s_ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
        :param seq_id_thr: Minimum sequence identity to be considered as match
        :type seq_id_thr: :class:`float`
        :param min_aln_length: Minimum number of aligned columns to be considered
                               as match. Only has an effect if *seq_id_thr* > 0!
        :type min_aln_length: :class:`int`    
        :returns: Tuple with two elements. 1) index of sequence in *ref_seqs* to
                  which *s* can be mapped 2) Pairwise sequence alignment with 
                  sequence from *ref_seqs* as first sequence. Both elements are
                  None if no mapping can be found or if thresholds are not 
                  fulfilled for any alignment.
        :raises: :class:`RuntimeError` if mapping is ambiguous, i.e. *s*
                 successfully maps to more than one sequence in *ref_seqs* 
        """
        scored_alns = list()
        for ref_idx, ref_seq in enumerate(ref_seqs):
            if ref_types[ref_idx] == s_type:
                if self.resnum_alignments:
                    aln = self.ResNumAlign(ref_seq, s, self.target, s_ent)
                else:
                    aln = self.NWAlign(ref_seq, s, s_type)
                seqid, n_tot, n_aligned = _GetAlnPropsTwo(aln)
                if seq_id_thr > 0:
                    if seqid >= seq_id_thr and n_aligned >= min_aln_length:
                        fraction_covered = float(n_aligned)/n_tot
                        score = seqid * fraction_covered
                        scored_alns.append((score, ref_idx, aln))
                else:
                    fraction_covered = float(n_aligned)/n_tot
                    score = seqid * fraction_covered
                    scored_alns.append((score, ref_idx, aln))
    
        if len(scored_alns) == 0:
            return (None, None) # no mapping possible...
    
        scored_alns = sorted(scored_alns, key=lambda x: x[0], reverse=True)
        return (scored_alns[0][1], scored_alns[0][2])

def _GetAlnPropsTwo(aln):
    """Returns basic properties of *aln* version two...

    :param aln: Alignment to compute properties
    :type aln: :class:`seq.AlignmentHandle`
    :returns: Tuple with 3 elements. 1) sequence identity in range [0, 100] 
              considering aligned columns 2) Number of non gap characters in
              first sequence 3) Number of aligned characters.
    """
    assert(aln.GetCount() == 2)
    n_tot = sum([1 for col in aln if col[0] != '-'])
    n_aligned = sum([1 for col in aln if (col[0] != '-' and col[1] != '-')])
    return (seq.alg.SequenceIdentity(aln), n_tot, n_aligned) 

def _GetAlnPropsOne(aln):
    
    """Returns basic properties of *aln* version one...

    :param aln: Alignment to compute properties
    :type aln: :class:`seq.AlignmentHandle`
    :returns: Tuple with 2 elements. 1) sequence identify in range [0, 100] 
              considering aligned columns 2) Number of aligned columns.
    """
    assert(aln.GetCount() == 2)
    seqid = seq.alg.SequenceIdentity(aln)
    n_aligned = sum([1 for col in aln if (col[0] != '-' and col[1] != '-')])
    return (seqid, n_aligned) 

def _GetRefMdlAlns(ref_chem_groups, ref_chem_group_msas, mdl_chem_groups,
                   mdl_chem_group_alns, pairs=None):
    """ Get all possible ref/mdl chain alignments given chem group mapping

    :param ref_chem_groups: :attr:`ChainMapper.chem_groups`
    :type ref_chem_groups: :class:`list` of :class:`list` of :class:`str`
    :param ref_chem_group_msas: :attr:`ChainMapper.chem_group_alignments`
    :type ref_chem_group_msas: :class:`ost.seq.AlignmentList`
    :param mdl_chem_groups: Groups of model chains that are mapped to
                            *ref_chem_groups*. Return value of
                            :func:`ChainMapper.GetChemMapping`.
    :type mdl_chem_groups: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chem_group_alns: A pairwise sequence alignment for every chain
                                in *mdl_chem_groups* that aligns these sequences
                                to the respective reference sequence.
                                Return values of
                                :func:`ChainMapper.GetChemMapping`.
    :type mdl_chem_group_alns: :class:`list` of :class:`ost.seq.AlignmentList`
    :param pairs: Pro param - restrict return dict to specified pairs. A set of
                  tuples in form (<trg_ch>, <mdl_ch>)
    :type pairs: :class:`set`
    :returns: A dictionary holding all possible ref/mdl chain alignments. Keys
              in that dictionary are tuples of the form (ref_ch, mdl_ch) and
              values are the respective pairwise alignments with first sequence
              being from ref, the second from mdl.
    """
    # alignment of each model chain to chem_group reference sequence
    mdl_alns = dict()
    for alns in mdl_chem_group_alns:
        for aln in alns:
            mdl_chain_name = aln.GetSequence(1).GetName()
            mdl_alns[mdl_chain_name] = aln

    # generate all alignments between ref/mdl chain atomseqs that we will
    # ever observe
    ref_mdl_alns = dict()
    for ref_chains, mdl_chains, ref_aln in zip(ref_chem_groups, mdl_chem_groups,
                                               ref_chem_group_msas):
        for ref_ch in ref_chains:
            for mdl_ch in mdl_chains:
                if pairs is not None and (ref_ch, mdl_ch) not in pairs:
                    continue
                # obtain alignments of mdl and ref chains towards chem
                # group ref sequence and merge them
                aln_list = seq.AlignmentList()
                
                # do ref aln
                ############
                # reference sequence
                s1 = ref_aln.GetSequence(0)
                # ATOMSEQ of ref_ch
                s2 = ref_aln.GetSequence(1+ref_chains.index(ref_ch))
                aln_list.append(seq.CreateAlignment(s1, s2))

                # do mdl aln
                ############
                aln_list.append(mdl_alns[mdl_ch])

                # merge
                #######
                ref_seq = seq.CreateSequence(s1.GetName(),
                                             s1.GetGaplessString())
                merged_aln = seq.alg.MergePairwiseAlignments(aln_list,
                                                             ref_seq)
                # merged_aln:
                # seq1: ref seq of chem group
                # seq2: seq of ref chain
                # seq3: seq of mdl chain
                # => we need the alignment between seq2 and seq3
                s2 = merged_aln.GetSequence(1)
                s3 = merged_aln.GetSequence(2)
                # cut leading and trailing gap columns
                a = 0 # number of leading gap columns
                for idx in range(len(s2)):
                    if s2[idx] != '-' or s3[idx] != '-':
                        break
                    a += 1
                b = 0 # number of trailing gap columns
                for idx in reversed(range(len(s2))):
                    if s2[idx] != '-' or s3[idx] != '-':
                        break
                    b += 1
                s2 = seq.CreateSequence(s2.GetName(), s2[a: len(s2)-b])
                s3 = seq.CreateSequence(s3.GetName(), s3[a: len(s3)-b])
                ref_mdl_alns[(ref_ch, mdl_ch)] = seq.CreateAlignment(s2, s3)

    return ref_mdl_alns

def _CheckOneToOneMapping(ref_chains, mdl_chains):
    """ Checks whether we already have a perfect one to one mapping

    That means each list in *ref_chains* has exactly one element and each
    list in *mdl_chains* has either one element (it's mapped) or is empty
    (ref chain has no mapped mdl chain). Returns None if no such mapping
    can be found.

    :param ref_chains: corresponds to :attr:`ChainMapper.chem_groups`
    :type ref_chains: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chains: mdl chains mapped to chem groups in *ref_chains*, i.e.
                       the return value of :func:`ChainMapper.GetChemMapping`
    :type mdl_chains: class:`list` of :class:`list` of :class:`str`
    :returns: A :class:`list` of :class:`list` if a one to one mapping is found,
              None otherwise
    """
    only_one_to_one = True
    one_to_one = list()
    for ref, mdl in zip(ref_chains, mdl_chains):
        if len(ref) == 1 and len(mdl) == 1:
            one_to_one.append(mdl)
        elif len(ref) == 1 and len(mdl) == 0:
            one_to_one.append([None])
        else:
            only_one_to_one = False
            break
    if only_one_to_one:
        return one_to_one
    else:
        return None

class _lDDTGreedySearcher(bb_lddt.BBlDDTScorer):
    def __init__(self, ref, mdl, ref_chem_groups, mdl_chem_groups,
                 ref_mdl_alns, inclusion_radius = 15.0,
                 thresholds = [0.5, 1.0, 2.0, 4.0],
                 steep_opt_rate = None):

        """ Greedy extension of already existing but incomplete chain mappings
        """
        super().__init__(ref, ref_chem_groups, mdl, ref_mdl_alns,
                         dist_thresh=inclusion_radius,
                         dist_diff_thresholds=thresholds)

        self.mdl_chem_groups = mdl_chem_groups
        self.steep_opt_rate = steep_opt_rate

        self.neighbors = {k: set() for k in self.trg.chain_names}
        for interface in self.trg.interacting_chains:
            self.neighbors[interface[0]].add(interface[1])
            self.neighbors[interface[1]].add(interface[0])

        assert(len(ref_chem_groups) == len(mdl_chem_groups))
        self.ref_chem_groups = ref_chem_groups
        self.mdl_chem_groups = mdl_chem_groups
        self.ref_ch_group_mapper = dict()
        self.mdl_ch_group_mapper = dict()
        for g_idx, (ref_g, mdl_g) in enumerate(zip(ref_chem_groups,
                                                   mdl_chem_groups)):
            for ch in ref_g:
                self.ref_ch_group_mapper[ch] = g_idx
            for ch in mdl_g:
                self.mdl_ch_group_mapper[ch] = g_idx

        # keep track of mdl chains that potentially give lDDT contributions,
        # i.e. they have locations within inclusion_radius + max(thresholds)
        self.mdl_neighbors = {k: set() for k in self.mdl.chain_names}
        for interface in self.mdl.potentially_contributing_chains:
            self.mdl_neighbors[interface[0]].add(interface[1])
            self.mdl_neighbors[interface[1]].add(interface[0])

    def ExtendMapping(self, mapping, max_ext = None):

        if len(mapping) == 0:
            raise RuntimError("Mapping must contain a starting point")

        # Ref chains onto which we can map. The algorithm starts with a mapping
        # on ref_ch. From there we can start to expand to connected neighbors.
        # All neighbors that we can reach from the already mapped chains are
        # stored in this set which will be updated during runtime
        map_targets = set()
        for ref_ch in mapping.keys():
            map_targets.update(self.neighbors[ref_ch])

        # remove the already mapped chains
        for ref_ch in mapping.keys():
            map_targets.discard(ref_ch)

        if len(map_targets) == 0:
            return mapping # nothing to extend

        # keep track of what model chains are not yet mapped for each chem group
        free_mdl_chains = list()
        for chem_group in self.mdl_chem_groups:
            tmp = [x for x in chem_group if x not in mapping.values()]
            free_mdl_chains.append(set(tmp))

        # keep track of what ref chains got a mapping
        newly_mapped_ref_chains = list()

        something_happened = True
        while something_happened:
            something_happened=False

            if self.steep_opt_rate is not None:
                n_chains = len(newly_mapped_ref_chains)
                if n_chains > 0 and n_chains % self.steep_opt_rate == 0:
                    mapping = self._SteepOpt(mapping, newly_mapped_ref_chains)

            if max_ext is not None and len(newly_mapped_ref_chains) >= max_ext:
                break

            max_n = 0
            max_mapping = None
            for ref_ch in map_targets:
                chem_group_idx = self.ref_ch_group_mapper[ref_ch]
                for mdl_ch in free_mdl_chains[chem_group_idx]:
                    # single chain score
                    n_single = self._NSCConserved(ref_ch, mdl_ch).sum()
                    # scores towards neighbors that are already mapped
                    n_inter = 0
                    for neighbor in self.neighbors[ref_ch]:
                        if neighbor in mapping and mapping[neighbor] in \
                        self.mdl_neighbors[mdl_ch]:
                            n_inter += self._NPairConserved((ref_ch, neighbor),
                                                            (mdl_ch, mapping[neighbor])).sum()
                    n = n_single + n_inter

                    if n_inter > 0 and n > max_n:
                        # Only accept a new solution if its actually connected
                        # i.e. n_inter > 0. Otherwise we could just map a big
                        # fat mdl chain sitting somewhere in Nirvana
                        max_n = n
                        max_mapping = (ref_ch, mdl_ch)
     
            if max_n > 0:
                something_happened = True
                # assign new found mapping
                mapping[max_mapping[0]] = max_mapping[1]

                # add all neighboring chains to map targets as they are now
                # reachable
                for neighbor in self.neighbors[max_mapping[0]]:
                    if neighbor not in mapping:
                        map_targets.add(neighbor)

                # remove the ref chain from map targets
                map_targets.remove(max_mapping[0])

                # remove the mdl chain from free_mdl_chains - its taken...
                chem_group_idx = self.ref_ch_group_mapper[max_mapping[0]]
                free_mdl_chains[chem_group_idx].remove(max_mapping[1])

                # keep track of what ref chains got a mapping
                newly_mapped_ref_chains.append(max_mapping[0])

        return mapping

    def _SteepOpt(self, mapping, chains_to_optimize=None):

        # just optimize ALL ref chains if nothing specified
        if chains_to_optimize is None:
            chains_to_optimize = mapping.keys()

        # make sure that we only have ref chains which are actually mapped
        ref_chains = [x for x in chains_to_optimize if mapping[x] is not None]

        # group ref chains to be optimized into chem groups
        tmp = dict()
        for ch in ref_chains:
            chem_group_idx = self.ref_ch_group_mapper[ch] 
            if chem_group_idx in tmp:
                tmp[chem_group_idx].append(ch)
            else:
                tmp[chem_group_idx] = [ch]
        chem_groups = list(tmp.values())

        # try all possible mapping swaps. Swaps that improve the score are
        # immediately accepted and we start all over again
        current_lddt = self.FromFlatMapping(mapping)
        something_happened = True
        while something_happened:
            something_happened = False
            for chem_group in chem_groups:
                if something_happened:
                    break
                for ch1, ch2 in itertools.combinations(chem_group, 2):
                    swapped_mapping = dict(mapping)
                    swapped_mapping[ch1] = mapping[ch2]
                    swapped_mapping[ch2] = mapping[ch1]
                    score = self.FromFlatMapping(swapped_mapping)
                    if score > current_lddt:
                        something_happened = True
                        mapping = swapped_mapping
                        current_lddt = score
                        break        
        return mapping


def _lDDTNaive(trg, mdl, inclusion_radius, thresholds, chem_groups,
               chem_mapping, ref_mdl_alns, n_max_naive):
    """ Naively iterates all possible chain mappings and returns the best
    """
    best_mapping = None
    best_lddt = -1.0

    # Setup scoring
    lddt_scorer = bb_lddt.BBlDDTScorer(trg, chem_groups, mdl, ref_mdl_alns,
                                       dist_thresh=inclusion_radius,
                                       dist_diff_thresholds=thresholds)
    for mapping in _ChainMappings(chem_groups, chem_mapping, n_max_naive):
        lDDT = lddt_scorer.Score(mapping, check=False)
        if lDDT > best_lddt:
            best_mapping = mapping
            best_lddt = lDDT

    return (best_mapping, best_lddt)


def _GetSeeds(ref_chem_groups, mdl_chem_groups, mapped_ref_chains = set(),
              mapped_mdl_chains = set()):
    seeds = list()
    for ref_chains, mdl_chains in zip(ref_chem_groups,
                                      mdl_chem_groups):
        for ref_ch in ref_chains:
            if ref_ch not in mapped_ref_chains:
                for mdl_ch in mdl_chains:
                    if mdl_ch not in mapped_mdl_chains:
                        seeds.append((ref_ch, mdl_ch))
    return seeds


def _lDDTGreedyFast(the_greed):

    something_happened = True
    mapping = dict()

    while something_happened:
        something_happened = False
        seeds = _GetSeeds(the_greed.ref_chem_groups,
                          the_greed.mdl_chem_groups,
                          mapped_ref_chains = set(mapping.keys()),
                          mapped_mdl_chains = set(mapping.values()))
        # search for best scoring starting point
        n_best = 0
        best_seed = None
        for seed in seeds:
            n = the_greed._NSCConserved(seed[0], seed[1]).sum()
            if n > n_best:
                n_best = n
                best_seed = seed
        if n_best == 0:
            break # no proper seed found anymore...
        # add seed to mapping and start the greed
        mapping[best_seed[0]] = best_seed[1]
        mapping = the_greed.ExtendMapping(mapping)
        something_happened = True

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


def _lDDTGreedyFull(the_greed):
    """ Uses each reference chain as starting point for expansion
    """

    seeds = _GetSeeds(the_greed.ref_chem_groups, the_greed.mdl_chem_groups)
    best_overall_score = -1.0
    best_overall_mapping = dict()

    for seed in seeds:

        # do initial extension
        mapping = the_greed.ExtendMapping({seed[0]: seed[1]})

        # repeat the process until we have a full mapping
        something_happened = True
        while something_happened:
            something_happened = False
            remnant_seeds = _GetSeeds(the_greed.ref_chem_groups,
                                      the_greed.mdl_chem_groups,
                                      mapped_ref_chains = set(mapping.keys()),
                                      mapped_mdl_chains = set(mapping.values()))
            if len(remnant_seeds) > 0:
                # still more mapping to be done
                best_score = -1.0
                best_mapping = None
                for remnant_seed in remnant_seeds:
                    tmp_mapping = dict(mapping)
                    tmp_mapping[remnant_seed[0]] = remnant_seed[1]
                    tmp_mapping = the_greed.ExtendMapping(tmp_mapping)
                    score = the_greed.FromFlatMapping(tmp_mapping)
                    if score > best_score:
                        best_score = score
                        best_mapping = tmp_mapping
                if best_mapping is not None:
                    something_happened = True
                    mapping = best_mapping

        score = the_greed.FromFlatMapping(mapping)
        if score > best_overall_score:
            best_overall_score = score
            best_overall_mapping = mapping

    mapping = best_overall_mapping

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


def _lDDTGreedyBlock(the_greed, seed_size, blocks_per_chem_group):
    """ try multiple seeds, i.e. try all ref/mdl chain combinations within the
    respective chem groups and compute single chain lDDTs. The
    *blocks_per_chem_group* best scoring ones are extend by *seed_size* chains
    and the best scoring one is exhaustively extended.
    """

    if seed_size is None or seed_size < 1:
        raise RuntimeError(f"seed_size must be an int >= 1 (got {seed_size})")

    if blocks_per_chem_group is None or blocks_per_chem_group < 1:
        raise RuntimeError(f"blocks_per_chem_group must be an int >= 1 "
                           f"(got {blocks_per_chem_group})")

    max_ext = seed_size - 1 #  -1 => start seed already has size 1

    ref_chem_groups = copy.deepcopy(the_greed.ref_chem_groups)
    mdl_chem_groups = copy.deepcopy(the_greed.mdl_chem_groups)

    mapping = dict()
    something_happened = True
    while something_happened:
        something_happened = False
        starting_blocks = list()
        for ref_chains, mdl_chains in zip(ref_chem_groups, mdl_chem_groups):
            if len(mdl_chains) == 0:
                continue # nothing to map
            ref_chains_copy = list(ref_chains)
            for i in range(blocks_per_chem_group):
                if len(ref_chains_copy) == 0:
                    break
                seeds = list()
                for ref_ch in ref_chains_copy:
                    seeds += [(ref_ch, mdl_ch) for mdl_ch in mdl_chains]
                # extend starting seeds to *seed_size* and retain best scoring
                # block for further extension
                best_score = -1.0
                best_mapping = None
                best_seed = None
                for s in seeds:
                    seed = dict(mapping)
                    seed.update({s[0]: s[1]})  
                    seed = the_greed.ExtendMapping(seed, max_ext = max_ext)
                    seed_lddt = the_greed.FromFlatMapping(seed)
                    if seed_lddt > best_score:
                        best_score = seed_lddt
                        best_mapping = seed
                        best_seed = s
                if best_mapping != None:
                    starting_blocks.append(best_mapping)
                    if best_seed[0] in ref_chains_copy:
                        # remove that ref chain to enforce diversity
                        ref_chains_copy.remove(best_seed[0])

        # fully expand initial starting blocks
        best_lddt = 0.0
        best_mapping = None
        for seed in starting_blocks:
            seed = the_greed.ExtendMapping(seed)
            seed_lddt = the_greed.FromFlatMapping(seed)
            if seed_lddt > best_lddt:
                best_lddt = seed_lddt
                best_mapping = seed

        if best_lddt == 0.0:
            break # no proper mapping found anymore

        something_happened = True
        mapping.update(best_mapping)
        for ref_ch, mdl_ch in best_mapping.items():
            for group_idx in range(len(ref_chem_groups)):
                if ref_ch in ref_chem_groups[group_idx]:
                    ref_chem_groups[group_idx].remove(ref_ch)
                if mdl_ch in mdl_chem_groups[group_idx]:
                    mdl_chem_groups[group_idx].remove(mdl_ch)

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


class _QSScoreGreedySearcher(qsscore.QSScorer):
    def __init__(self, ref, mdl, ref_chem_groups, mdl_chem_groups,
                 ref_mdl_alns, contact_d = 12.0,
                 steep_opt_rate = None, greedy_prune_contact_map=False):
        """ Greedy extension of already existing but incomplete chain mappings
        """
        super().__init__(ref, ref_chem_groups, mdl, ref_mdl_alns,
                         contact_d = contact_d)
        self.ref = ref
        self.mdl = mdl
        self.ref_mdl_alns = ref_mdl_alns
        self.steep_opt_rate = steep_opt_rate

        if greedy_prune_contact_map:
            self.neighbors = {k: set() for k in self.qsent1.chain_names}
            for p in self.qsent1.interacting_chains:
                if np.count_nonzero(self.qsent1.PairDist(p[0], p[1])<=8) >= 3:
                    self.neighbors[p[0]].add(p[1])
                    self.neighbors[p[1]].add(p[0])

            self.mdl_neighbors = {k: set() for k in self.qsent2.chain_names}
            for p in self.qsent2.interacting_chains:
                if np.count_nonzero(self.qsent2.PairDist(p[0], p[1])<=8) >= 3:
                    self.mdl_neighbors[p[0]].add(p[1])
                    self.mdl_neighbors[p[1]].add(p[0])


        else:
            self.neighbors = {k: set() for k in self.qsent1.chain_names}
            for p in self.qsent1.interacting_chains:
                self.neighbors[p[0]].add(p[1])
                self.neighbors[p[1]].add(p[0])

            self.mdl_neighbors = {k: set() for k in self.qsent2.chain_names}
            for p in self.qsent2.interacting_chains:
                self.mdl_neighbors[p[0]].add(p[1])
                self.mdl_neighbors[p[1]].add(p[0])

        assert(len(ref_chem_groups) == len(mdl_chem_groups))
        self.ref_chem_groups = ref_chem_groups
        self.mdl_chem_groups = mdl_chem_groups
        self.ref_ch_group_mapper = dict()
        self.mdl_ch_group_mapper = dict()
        for g_idx, (ref_g, mdl_g) in enumerate(zip(ref_chem_groups,
                                                   mdl_chem_groups)):
            for ch in ref_g:
                self.ref_ch_group_mapper[ch] = g_idx
            for ch in mdl_g:
                self.mdl_ch_group_mapper[ch] = g_idx

        # cache for lDDT based single chain conserved contacts
        # used to identify starting points for further extension by QS score
        # key: tuple (ref_ch, mdl_ch) value: number of conserved contacts
        self.single_chain_scorer = dict()
        self.single_chain_cache = dict()
        for ch in self.ref.chains:
            single_chain_ref = _CSel(self.ref, [ch.GetName()])
            self.single_chain_scorer[ch.GetName()] = \
            lddt.lDDTScorer(single_chain_ref, bb_only = True)

    def SCCounts(self, ref_ch, mdl_ch):
        if not (ref_ch, mdl_ch) in self.single_chain_cache:
            alns = dict()
            alns[mdl_ch] = self.ref_mdl_alns[(ref_ch, mdl_ch)]
            mdl_sel = _CSel(self.mdl, [mdl_ch])
            s = self.single_chain_scorer[ref_ch]
            _,_,_,conserved,_,_,_ = s.lDDT(mdl_sel,
                                           residue_mapping=alns,
                                           return_dist_test=True,
                                           no_interchain=True,
                                           chain_mapping={mdl_ch: ref_ch},
                                           check_resnames=False)
            self.single_chain_cache[(ref_ch, mdl_ch)] = conserved
        return self.single_chain_cache[(ref_ch, mdl_ch)]

    def ExtendMapping(self, mapping, max_ext = None):

        if len(mapping) == 0:
            raise RuntimError("Mapping must contain a starting point")

        # Ref chains onto which we can map. The algorithm starts with a mapping
        # on ref_ch. From there we can start to expand to connected neighbors.
        # All neighbors that we can reach from the already mapped chains are
        # stored in this set which will be updated during runtime
        map_targets = set()
        for ref_ch in mapping.keys():
            map_targets.update(self.neighbors[ref_ch])

        # remove the already mapped chains
        for ref_ch in mapping.keys():
            map_targets.discard(ref_ch)

        if len(map_targets) == 0:
            return mapping # nothing to extend

        # keep track of what model chains are not yet mapped for each chem group
        free_mdl_chains = list()
        for chem_group in self.mdl_chem_groups:
            tmp = [x for x in chem_group if x not in mapping.values()]
            free_mdl_chains.append(set(tmp))

        # keep track of what ref chains got a mapping
        newly_mapped_ref_chains = list()

        something_happened = True
        while something_happened:
            something_happened=False

            if self.steep_opt_rate is not None:
                n_chains = len(newly_mapped_ref_chains)
                if n_chains > 0 and n_chains % self.steep_opt_rate == 0:
                    mapping = self._SteepOpt(mapping, newly_mapped_ref_chains)

            if max_ext is not None and len(newly_mapped_ref_chains) >= max_ext:
                break

            score_result = self.FromFlatMapping(mapping)
            old_score = score_result.QS_global
            nominator = score_result.weighted_scores
            denominator = score_result.weight_sum + score_result.weight_extra_all

            max_diff = 0.0
            max_mapping = None
            for ref_ch in map_targets:
                chem_group_idx = self.ref_ch_group_mapper[ref_ch]
                for mdl_ch in free_mdl_chains[chem_group_idx]:
                    # we're not computing full QS-score here, we directly hack
                    # into the QS-score formula to compute a diff
                    nominator_diff = 0.0
                    denominator_diff = 0.0
                    for neighbor in self.neighbors[ref_ch]:
                        if neighbor in mapping and mapping[neighbor] in \
                        self.mdl_neighbors[mdl_ch]:
                            # it's a newly added interface if (ref_ch, mdl_ch)
                            # are added to mapping
                            int1 = (ref_ch, neighbor)
                            int2 = (mdl_ch, mapping[neighbor])
                            a, b, c, d = self._MappedInterfaceScores(int1, int2)
                            nominator_diff += a # weighted_scores
                            denominator_diff += b # weight_sum
                            denominator_diff += d # weight_extra_all
                            # the respective interface penalties are subtracted
                            # from denominator
                            denominator_diff -= self._InterfacePenalty1(int1)
                            denominator_diff -= self._InterfacePenalty2(int2)

                    if nominator_diff > 0:
                        # Only accept a new solution if its actually connected
                        # i.e. nominator_diff > 0.
                        new_nominator = nominator + nominator_diff
                        new_denominator = denominator + denominator_diff
                        new_score = 0.0
                        if new_denominator != 0.0:
                            new_score = new_nominator/new_denominator
                        diff = new_score - old_score
                        if diff > max_diff:
                            max_diff = diff
                            max_mapping = (ref_ch, mdl_ch)
     
            if max_mapping is not None:
                something_happened = True
                # assign new found mapping
                mapping[max_mapping[0]] = max_mapping[1]

                # add all neighboring chains to map targets as they are now
                # reachable
                for neighbor in self.neighbors[max_mapping[0]]:
                    if neighbor not in mapping:
                        map_targets.add(neighbor)

                # remove the ref chain from map targets
                map_targets.remove(max_mapping[0])

                # remove the mdl chain from free_mdl_chains - its taken...
                chem_group_idx = self.ref_ch_group_mapper[max_mapping[0]]
                free_mdl_chains[chem_group_idx].remove(max_mapping[1])

                # keep track of what ref chains got a mapping
                newly_mapped_ref_chains.append(max_mapping[0])

        return mapping

    def _SteepOpt(self, mapping, chains_to_optimize=None):

        # just optimize ALL ref chains if nothing specified
        if chains_to_optimize is None:
            chains_to_optimize = mapping.keys()

        # make sure that we only have ref chains which are actually mapped
        ref_chains = [x for x in chains_to_optimize if mapping[x] is not None]

        # group ref chains to be optimized into chem groups
        tmp = dict()
        for ch in ref_chains:
            chem_group_idx = self.ref_ch_group_mapper[ch] 
            if chem_group_idx in tmp:
                tmp[chem_group_idx].append(ch)
            else:
                tmp[chem_group_idx] = [ch]
        chem_groups = list(tmp.values())

        # try all possible mapping swaps. Swaps that improve the score are
        # immediately accepted and we start all over again
        score_result = self.FromFlatMapping(mapping)
        current_score = score_result.QS_global
        something_happened = True
        while something_happened:
            something_happened = False
            for chem_group in chem_groups:
                if something_happened:
                    break
                for ch1, ch2 in itertools.combinations(chem_group, 2):
                    swapped_mapping = dict(mapping)
                    swapped_mapping[ch1] = mapping[ch2]
                    swapped_mapping[ch2] = mapping[ch1]
                    score_result = self.FromFlatMapping(swapped_mapping)
                    if score_result.QS_global > current_score:
                        something_happened = True
                        mapping = swapped_mapping
                        current_score = score_result.QS_global
                        break        
        return mapping


def _QSScoreNaive(trg, mdl, chem_groups, chem_mapping, ref_mdl_alns, contact_d,
                  n_max_naive):
    best_mapping = None
    best_score = -1.0
    # qs_scorer implements caching, score calculation is thus as fast as it gets
    # you'll just hit a wall when the number of possible mappings becomes large
    qs_scorer = qsscore.QSScorer(trg, chem_groups, mdl, ref_mdl_alns)
    for mapping in _ChainMappings(chem_groups, chem_mapping, n_max_naive):
        score_result = qs_scorer.Score(mapping, check=False)
        if score_result.QS_global > best_score:
            best_mapping = mapping
            best_score = score_result.QS_global
    return (best_mapping, best_score)


def _QSScoreGreedyFast(the_greed):

    something_happened = True
    mapping = dict()
    while something_happened:
        something_happened = False
        # search for best scoring starting point, we're using lDDT here
        n_best = 0
        best_seed = None
        seeds = _GetSeeds(the_greed.ref_chem_groups,
                          the_greed.mdl_chem_groups,
                          mapped_ref_chains = set(mapping.keys()),
                          mapped_mdl_chains = set(mapping.values()))
        for seed in seeds:
            n = the_greed.SCCounts(seed[0], seed[1])
            if n > n_best:
                n_best = n
                best_seed = seed
        if n_best == 0:
            break # no proper seed found anymore...
        # add seed to mapping and start the greed
        mapping[best_seed[0]] = best_seed[1]
        mapping = the_greed.ExtendMapping(mapping)
        something_happened = True

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


def _QSScoreGreedyFull(the_greed):
    """ Uses each reference chain as starting point for expansion
    """

    seeds = _GetSeeds(the_greed.ref_chem_groups, the_greed.mdl_chem_groups)
    best_overall_score = -1.0
    best_overall_mapping = dict()

    for seed in seeds:

        # do initial extension
        mapping = the_greed.ExtendMapping({seed[0]: seed[1]})

        # repeat the process until we have a full mapping
        something_happened = True
        while something_happened:
            something_happened = False
            remnant_seeds = _GetSeeds(the_greed.ref_chem_groups,
                                      the_greed.mdl_chem_groups,
                                      mapped_ref_chains = set(mapping.keys()),
                                      mapped_mdl_chains = set(mapping.values()))
            if len(remnant_seeds) > 0:
                # still more mapping to be done
                best_score = -1.0
                best_mapping = None
                for remnant_seed in remnant_seeds:
                    tmp_mapping = dict(mapping)
                    tmp_mapping[remnant_seed[0]] = remnant_seed[1]
                    tmp_mapping = the_greed.ExtendMapping(tmp_mapping)
                    score_result = the_greed.FromFlatMapping(tmp_mapping)
                    if score_result.QS_global > best_score:
                        best_score = score_result.QS_global
                        best_mapping = tmp_mapping
                if best_mapping is not None:
                    something_happened = True
                    mapping = best_mapping

        score_result = the_greed.FromFlatMapping(mapping)
        if score_result.QS_global > best_overall_score:
            best_overall_score = score_result.QS_global
            best_overall_mapping = mapping

    mapping = best_overall_mapping

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping


def _QSScoreGreedyBlock(the_greed, seed_size, blocks_per_chem_group):
    """ try multiple seeds, i.e. try all ref/mdl chain combinations within the
    respective chem groups and compute single chain lDDTs. The
    *blocks_per_chem_group* best scoring ones are extend by *seed_size* chains
    and the best scoring one with respect to QS score is exhaustively extended.
    """

    if seed_size is None or seed_size < 1:
        raise RuntimeError(f"seed_size must be an int >= 1 (got {seed_size})")

    if blocks_per_chem_group is None or blocks_per_chem_group < 1:
        raise RuntimeError(f"blocks_per_chem_group must be an int >= 1 "
                           f"(got {blocks_per_chem_group})")

    max_ext = seed_size - 1 #  -1 => start seed already has size 1

    ref_chem_groups = copy.deepcopy(the_greed.ref_chem_groups)
    mdl_chem_groups = copy.deepcopy(the_greed.mdl_chem_groups)

    mapping = dict()

    something_happened = True
    while something_happened:
        something_happened = False
        starting_blocks = list()
        for ref_chains, mdl_chains in zip(ref_chem_groups, mdl_chem_groups):
            if len(mdl_chains) == 0:
                continue # nothing to map
            ref_chains_copy = list(ref_chains)
            for i in range(blocks_per_chem_group):
                if len(ref_chains_copy) == 0:
                    break
                seeds = list()
                for ref_ch in ref_chains_copy:
                    seeds += [(ref_ch, mdl_ch) for mdl_ch in mdl_chains]
                # extend starting seeds to *seed_size* and retain best scoring block
                # for further extension
                best_score = -1.0
                best_mapping = None
                best_seed = None
                for s in seeds:
                    seed = dict(mapping)
                    seed.update({s[0]: s[1]})  
                    seed = the_greed.ExtendMapping(seed, max_ext = max_ext)
                    score_result = the_greed.FromFlatMapping(seed)
                    if score_result.QS_global > best_score:
                        best_score = score_result.QS_global
                        best_mapping = seed
                        best_seed = s
                if best_mapping != None:
                    starting_blocks.append(best_mapping)
                    if best_seed[0] in ref_chains_copy:
                        # remove selected ref chain to enforce diversity
                        ref_chains_copy.remove(best_seed[0])

        # fully expand initial starting blocks
        best_score = -1.0
        best_mapping = None
        for seed in starting_blocks:
            seed = the_greed.ExtendMapping(seed)
            score_result = the_greed.FromFlatMapping(seed)
            if score_result.QS_global > best_score:
                best_score = score_result.QS_global
                best_mapping = seed

        if best_mapping is not None and len(best_mapping) > len(mapping):
            # this even accepts extensions that lead to no increase in QS-score
            # at least they make sense from an lDDT perspective
            something_happened = True
            mapping.update(best_mapping)
            for ref_ch, mdl_ch in best_mapping.items():
                for group_idx in range(len(ref_chem_groups)):
                    if ref_ch in ref_chem_groups[group_idx]:
                        ref_chem_groups[group_idx].remove(ref_ch)
                    if mdl_ch in mdl_chem_groups[group_idx]:
                        mdl_chem_groups[group_idx].remove(mdl_ch)

    # translate mapping format and return
    final_mapping = list()
    for ref_chains in the_greed.ref_chem_groups:
        mapped_mdl_chains = list()
        for ref_ch in ref_chains:
            if ref_ch in mapping:
                mapped_mdl_chains.append(mapping[ref_ch])
            else:
                mapped_mdl_chains.append(None)
        final_mapping.append(mapped_mdl_chains)

    return final_mapping

def _SingleRigidRMSD(initial_transforms, initial_mappings, chem_groups,
                     chem_mapping, trg_group_pos, mdl_group_pos):
    """
    Takes initial transforms and sequentially adds chain pairs with lowest RMSD.
    The mapping from the transform that leads to lowest overall RMSD is
    returned.
    """
    best_mapping = dict()
    best_ssd = float("inf") # we're actually going for summed squared distances
                            # Since all positions have same lengths and we do a
                            # full mapping, lowest SSD has a guarantee of also
                            # being lowest RMSD
    for transform in initial_transforms:
        mapping = dict()
        mapped_mdl_chains = set()
        ssd = 0.0
        for trg_chains, mdl_chains, trg_pos, mdl_pos, in zip(chem_groups,
                                                             chem_mapping,
                                                             trg_group_pos,
                                                             mdl_group_pos):
            if len(trg_pos) == 0 or len(mdl_pos) == 0:
                continue # cannot compute valid rmsd
            ssds = list()
            t_mdl_pos = list()
            for m_pos in mdl_pos:
                t_m_pos = geom.Vec3List(m_pos)
                t_m_pos.ApplyTransform(transform)
                t_mdl_pos.append(t_m_pos)
            for t_pos, t in zip(trg_pos, trg_chains):
                for t_m_pos, m in zip(t_mdl_pos, mdl_chains):
                    ssd = t_pos.GetSummedSquaredDistances(t_m_pos)
                    ssds.append((ssd, (t,m)))
            ssds.sort()
            for item in ssds:
                p = item[1]
                if p[0] not in mapping and p[1] not in mapped_mdl_chains:
                    mapping[p[0]] = p[1]
                    mapped_mdl_chains.add(p[1])
                    ssd += item[0]

        if ssd < best_ssd:
            best_ssd = ssd
            best_mapping = mapping

    return best_mapping

def _IterativeRigidRMSD(initial_transforms, initial_mappings, chem_groups,
                        chem_mapping, trg_group_pos, mdl_group_pos):
    """ Takes initial transforms and sequentially adds chain pairs with
    lowest RMSD. With each added chain pair, the transform gets updated.
    Thus the naming iterative. The mapping from the initial transform that
    leads to best overall RMSD score is returned.
    """

    # to directly retrieve positions using chain names
    trg_pos_dict = dict()
    for trg_pos, trg_chains in zip(trg_group_pos, chem_groups):
        for t_pos, t in zip(trg_pos, trg_chains):
            trg_pos_dict[t] = t_pos
    mdl_pos_dict = dict()
    for mdl_pos, mdl_chains in zip(mdl_group_pos, chem_mapping):
        for m_pos, m in zip(mdl_pos, mdl_chains):
            mdl_pos_dict[m] = m_pos
        
    best_mapping = dict()
    best_rmsd = float("inf")
    for initial_transform, initial_mapping in zip(initial_transforms,
                                                  initial_mappings):
        mapping = {initial_mapping[0]: initial_mapping[1]}
        transform = geom.Mat4(initial_transform)
        mapped_trg_pos = geom.Vec3List(trg_pos_dict[initial_mapping[0]])
        mapped_mdl_pos = geom.Vec3List(mdl_pos_dict[initial_mapping[1]])

        # the following variables contain the chains which are
        # available for mapping
        trg_chain_groups = [set(group) for group in chem_groups]
        mdl_chain_groups = [set(group) for group in chem_mapping]

        # search and kick out inital mapping
        for group in trg_chain_groups:
            if initial_mapping[0] in group:
                group.remove(initial_mapping[0])
                break
        for group in mdl_chain_groups:
            if initial_mapping[1] in group:
                group.remove(initial_mapping[1])
                break

        something_happened = True
        while something_happened:
            # search for best mapping given current transform
            something_happened=False
            best_sc_mapping = None
            best_sc_group_idx = None
            best_sc_rmsd = float("inf")
            group_idx = 0
            for trg_chains, mdl_chains in zip(trg_chain_groups, mdl_chain_groups):
                for t in trg_chains:
                    t_pos = trg_pos_dict[t]
                    for m in mdl_chains:
                        m_pos = mdl_pos_dict[m]
                        t_m_pos = geom.Vec3List(m_pos)
                        t_m_pos.ApplyTransform(transform)
                        rmsd = t_pos.GetRMSD(t_m_pos)
                        if rmsd < best_sc_rmsd:
                            best_sc_rmsd = rmsd
                            best_sc_mapping = (t,m)
                            best_sc_group_idx = group_idx
                group_idx += 1

            if best_sc_mapping is not None:
                something_happened = True
                mapping[best_sc_mapping[0]] = best_sc_mapping[1]
                mapped_trg_pos.extend(trg_pos_dict[best_sc_mapping[0]])
                mapped_mdl_pos.extend(mdl_pos_dict[best_sc_mapping[1]])
                trg_chain_groups[best_sc_group_idx].remove(best_sc_mapping[0])
                mdl_chain_groups[best_sc_group_idx].remove(best_sc_mapping[1])

                transform = _GetSuperposition(mapped_mdl_pos, mapped_trg_pos,
                                              False).transformation

        # compute overall RMSD for current transform
        mapped_mdl_pos.ApplyTransform(transform)
        rmsd = mapped_trg_pos.GetRMSD(mapped_mdl_pos)

        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_mapping = mapping

    return best_mapping

def _NaiveRMSD(chem_groups, chem_mapping, trg_group_pos, mdl_group_pos,
               n_max_naive):

    # to directly retrieve positions using chain names
    trg_pos_dict = dict()
    for trg_pos, trg_chains in zip(trg_group_pos, chem_groups):
        for t_pos, t in zip(trg_pos, trg_chains):
            trg_pos_dict[t] = t_pos
    mdl_pos_dict = dict()
    for mdl_pos, mdl_chains in zip(mdl_group_pos, chem_mapping):
        for m_pos, m in zip(mdl_pos, mdl_chains):
            mdl_pos_dict[m] = m_pos
        
    best_mapping = dict()
    best_rmsd = float("inf")

    for mapping in _ChainMappings(chem_groups, chem_mapping, n_max_naive):
        trg_pos = geom.Vec3List()
        mdl_pos = geom.Vec3List()
        for trg_group, mdl_group in zip(chem_groups, mapping):
            for trg_ch, mdl_ch in zip(trg_group, mdl_group):
                if trg_ch is not None and mdl_ch is not None:
                    trg_pos.extend(trg_pos_dict[trg_ch])
                    mdl_pos.extend(mdl_pos_dict[mdl_ch])
        superpose_res = mol.alg.SuperposeSVD(mdl_pos, trg_pos)
        if superpose_res.rmsd < best_rmsd:
            best_rmsd = superpose_res.rmsd
            best_mapping = mapping

    # this is stupid...
    tmp = dict()
    for chem_group, mapping in zip(chem_groups, best_mapping):
        for trg_ch, mdl_ch in zip(chem_group, mapping):
            tmp[trg_ch] = mdl_ch

    return tmp


def _GetRefPos(trg, mdl, trg_msas, mdl_alns, max_pos = None):
    """ Extracts reference positions which are present in trg and mdl
    """

    # mdl_alns are pairwise, let's construct MSAs
    mdl_msas = list()
    for aln_list in mdl_alns:
        if len(aln_list) > 0:
            tmp = aln_list[0].GetSequence(0)
            ref_seq = seq.CreateSequence(tmp.GetName(), tmp.GetGaplessString())
            mdl_msas.append(seq.alg.MergePairwiseAlignments(aln_list, ref_seq))
        else:
            mdl_msas.append(seq.CreateAlignment())

    # fetch positions of all backbone atoms
    bb_trg = _GetBBPos(trg)
    bb_mdl = _GetBBPos(mdl)

    # there are the actual return values
    trg_pos = list()
    mdl_pos = list()

    for trg_msa, mdl_msa in zip(trg_msas, mdl_msas):

        if mdl_msa.GetCount() > 0:
            # make sure they have the same ref sequence (should be a given...)
            assert(trg_msa.GetSequence(0).GetGaplessString() == \
                   mdl_msa.GetSequence(0).GetGaplessString())
        else:
            # if mdl_msa is empty, i.e. no model chain maps to the chem group
            # represented by trg_msa, we just continue. The result will be
            # empty position lists added to trg_pos and mdl_pos.
            pass 

        # check which columns in MSAs are fully covered (indices relative to
        # first sequence)
        trg_indices = _GetFullyCoveredIndices(trg_msa)
        mdl_indices = _GetFullyCoveredIndices(mdl_msa)

        # get indices where both, mdl and trg, are fully covered
        indices = sorted(list(trg_indices.intersection(mdl_indices)))

        # subsample if necessary
        if max_pos is not None and len(indices) > max_pos:
            step = int(len(indices)/max_pos)
            indices = [indices[i] for i in range(0, len(indices), step)]

        # translate to column indices in the respective MSAs
        trg_indices = _RefIndicesToColumnIndices(trg_msa, indices)
        mdl_indices = _RefIndicesToColumnIndices(mdl_msa, indices)

        # extract positions
        trg_pos.append(list())
        mdl_pos.append(list())
        # first seq in trg_msa is ref sequence and does not refer to any
        # ATOMSEQ
        for s_idx in range(1, trg_msa.GetCount()):
            trg_pos[-1].append(_ExtractMSAPos(trg_msa, s_idx, trg_indices,
                                              bb_trg))
        # first seq in mdl_msa is ref sequence in trg and does not belong to mdl
        for s_idx in range(1, mdl_msa.GetCount()):
            mdl_pos[-1].append(_ExtractMSAPos(mdl_msa, s_idx, mdl_indices,
                                              bb_mdl))

    return (trg_pos, mdl_pos)

def _GetBBPos(ent):
    """ Helper for _GetRefPos

    Returns a dict with key: chain name and value: a geom.Vec3List of backbone
    atom positions. Assumes that either CA or C3' is always there.
    """
    bb = dict()
    for ch in ent.chains:
        l = geom.Vec3List()
        for r in ch.residues:
            a = r.FindAtom("CA")
            if not a.IsValid():
                a = r.FindAtom("C3'")
            l.append(a.GetPos())
        bb[ch.name] = l
    return bb

def _GetFullyCoveredIndices(msa):
    """ Helper for _GetRefPos

    Returns a set containing the indices relative to first sequence in msa which
    are fully covered in all other sequences

    --AA-A-A
    -BBBB-BB
    CCCC-C-C

    => (0,1,3)
    """
    indices = set()
    ref_idx = 0
    for col in msa:
        if sum([1 for olc in col if olc != '-']) == col.GetRowCount():
            indices.add(ref_idx)
        if col[0] != '-':
            ref_idx += 1
    return indices

def _RefIndicesToColumnIndices(msa, indices):
    """ Helper for _GetRefPos

    Returns a list of mapped indices. indices refer to non-gap one letter
    codes in the first msa sequence. The returnes mapped indices are translated
    to the according msa column indices
    """
    ref_idx = 0
    mapping = dict()
    for col_idx, col in enumerate(msa):
        if col[0] != '-':
            mapping[ref_idx] = col_idx
            ref_idx += 1
    return [mapping[i] for i in indices]

def _ExtractMSAPos(msa, s_idx, indices, bb):
    """ Helper for _GetRefPos

    Returns a geom.Vec3List containing positions refering to given msa sequence.
    => Chain with corresponding name is mapped onto sequence and the position of
    the first atom of each residue specified in indices is extracted.
    Indices refers to column indices in msa!
    """
    s = msa.GetSequence(s_idx)
    ch_bb = bb[s.GetName()]

    # sanity check
    assert(len(s.GetGaplessString()) == len(ch_bb))

    residue_idx = [s.GetResidueIndex(i) for i in indices]
    return geom.Vec3List([ch_bb[i] for i in residue_idx])


def _NChemGroupMappings(ref_chains, mdl_chains):
    """ Number of mappings within one chem group

    :param ref_chains: Reference chains
    :type ref_chains: :class:`list` of :class:`str`
    :param mdl_chains: Model chains that are mapped onto *ref_chains*
    :type mdl_chains: :class:`list` of :class:`str`
    :returns: Number of possible mappings of *mdl_chains* onto *ref_chains*
    """
    n_ref = len(ref_chains)
    n_mdl = len(mdl_chains)
    if n_ref == n_mdl:
        return factorial(n_ref)
    elif n_ref > n_mdl:
        n_choose_k = binom(n_ref, n_mdl)
        return n_choose_k * factorial(n_mdl)
    else:
        n_choose_k = binom(n_mdl, n_ref)
        return n_choose_k * factorial(n_ref)

def _NMappings(ref_chains, mdl_chains):
    """ Number of mappings for a full chem mapping

    :param ref_chains: Chem groups of reference
    :type ref_chains: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chains: Model chains that map onto those chem groups
    :type mdl_chains: :class:`list` of :class:`list` of :class:`str`
    :returns: Number of possible mappings of *mdl_chains* onto *ref_chains*
    """
    assert(len(ref_chains) == len(mdl_chains))
    n = 1
    for a,b in zip(ref_chains, mdl_chains):
        n *= _NChemGroupMappings(a,b)
    return n

def _NMappingsWithin(ref_chains, mdl_chains, max_mappings):
    """ Check whether total number of mappings is smaller than given maximum

    In principle the same as :func:`_NMappings` but it stops as soon as the
    maximum is hit.

    :param ref_chains: Chem groups of reference
    :type ref_chains: :class:`list` of :class:`list` of :class:`str`
    :param mdl_chains: Model chains that map onto those chem groups
    :type mdl_chains: :class:`list` of :class:`list` of :class:`str`
    :param max_mappings: Number of max allowed mappings
    :returns: Whether number of possible mappings of *mdl_chains* onto
              *ref_chains* is below or equal *max_mappings*.
    """
    assert(len(ref_chains) == len(mdl_chains))
    n = 1
    for a,b in zip(ref_chains, mdl_chains):
        n *= _NChemGroupMappings(a,b)
        if n > max_mappings:
            return False
    return True

def _RefSmallerGenerator(ref_chains, mdl_chains):
    """ Returns all possible ways to map mdl_chains onto ref_chains

    Specific for the case where len(ref_chains) < len(mdl_chains)
    """
    for c in itertools.combinations(mdl_chains, len(ref_chains)):
        for p in itertools.permutations(c):
            yield list(p)

def _RefLargerGenerator(ref_chains, mdl_chains):
    """ Returns all possible ways to map mdl_chains onto ref_chains

    Specific for the case where len(ref_chains) > len(mdl_chains)
    Ref chains without mapped mdl chain are assigned None
    """
    n_ref = len(ref_chains)
    n_mdl = len(mdl_chains)
    for c in itertools.combinations(range(n_ref), n_mdl):
        for p in itertools.permutations(mdl_chains):
            ret_list = [None] * n_ref
            for idx, ch in zip(c, p):
                ret_list[idx] = ch
            yield ret_list

def _RefEqualGenerator(ref_chains, mdl_chains):
    """ Returns all possible ways to map mdl_chains onto ref_chains

    Specific for the case where len(ref_chains) == len(mdl_chains)
    """
    for p in itertools.permutations(mdl_chains):
        yield list(p)

def _RefEmptyGenerator(ref_chains, mdl_chains):
    yield list()

def _ConcatIterators(iterators):
    for item in itertools.product(*iterators):
        yield list(item)

def _ChainMappings(ref_chains, mdl_chains, n_max=None):
    """Returns all possible ways to map *mdl_chains* onto fixed *ref_chains*

    :param ref_chains: List of list of chemically equivalent chains in reference
    :type ref_chains: :class:`list` of :class:`list`
    :param mdl_chains: Equally long list of list of chemically equivalent chains
                       in model that map on those ref chains.
    :type mdl_chains: :class:`list` of :class:`list`
    :param n_max: Aborts and raises :class:`RuntimeError` if max number of
                  mappings is above this threshold.
    :type n_max: :class:`int`
    :returns: Iterator over all possible mappings of *mdl_chains* onto fixed
              *ref_chains*. Potentially contains None as padding when number of
              model chains for a certain mapping is smaller than the according
              reference chains.
              Example: _ChainMappings([['A', 'B', 'C'], ['D', 'E']],
                                      [['x', 'y'], ['i', 'j']])
              gives an iterator over: [[['x', 'y', None], ['i', 'j']],
                                       [['x', 'y', None], ['j', 'i']],
                                       [['y', 'x', None], ['i', 'j']],
                                       [['y', 'x', None], ['j', 'i']],
                                       [['x', None, 'y'], ['i', 'j']],
                                       [['x', None, 'y'], ['j', 'i']],
                                       [['y', None, 'x'], ['i', 'j']],
                                       [['y', None, 'x'], ['j', 'i']],
                                       [[None, 'x', 'y'], ['i', 'j']],
                                       [[None, 'x', 'y'], ['j', 'i']],
                                       [[None, 'y', 'x'], ['i', 'j']],
                                       [[None, 'y', 'x'], ['j', 'i']]]
    """
    assert(len(ref_chains) == len(mdl_chains))

    if n_max is not None:
        if not _NMappingsWithin(ref_chains, mdl_chains, n_max):
            raise RuntimeError(f"Too many mappings. Max allowed: {n_max}")

    # one iterator per mapping representing all mdl combinations relative to
    # reference
    iterators = list()
    for ref, mdl in zip(ref_chains, mdl_chains):
        if len(ref) == 0:
            iterators.append(_RefEmptyGenerator(ref, mdl))
        elif len(ref) == len(mdl):
            iterators.append(_RefEqualGenerator(ref, mdl))
        elif len(ref) < len(mdl):
            iterators.append(_RefSmallerGenerator(ref, mdl))
        else:
            iterators.append(_RefLargerGenerator(ref, mdl))

    return _ConcatIterators(iterators)


def _GetSuperposition(pos_one, pos_two, iterative):
    """ Computes minimal RMSD superposition for pos_one onto pos_two

    :param pos_one: Positions that should be superposed onto *pos_two*
    :type pos_one: :class:`geom.Vec3List`
    :param pos_two: Reference positions
    :type pos_two: :class:`geom.Vec3List`
    :iterative: Whether iterative superposition should be used. Iterative
                potentially raises, uses standard superposition as fallback.
    :type iterative: :class:`bool`
    :returns: Transformation matrix to superpose *pos_one* onto *pos_two*
    :rtype: :class:`ost.mol.alg.SuperpositionResult`
    """
    res = None
    if iterative:
        try:
            res = mol.alg.IterativeSuperposeSVD(pos_one, pos_two)
        except:
            pass # triggers fallback below
    if res is None:
        res = mol.alg.SuperposeSVD(pos_one, pos_two)
    return res

# specify public interface
__all__ = ('ChainMapper', 'ReprResult', 'MappingResult')
