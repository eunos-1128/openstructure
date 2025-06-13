from contextlib import contextmanager
import numpy as np
import networkx

import ost
from ost import mol
from ost import conop
from ost import LogWarning, LogScript, LogInfo, LogVerbose, LogDebug, GetVerbosityLevel, PushVerbosityLevel, PopVerbosityLevel
from ost.mol.alg import chain_mapping

@contextmanager
def _SinkVerbosityLevel(n=1):
    """ Context manager to temporarily reduce the verbosity level by n.

    Example usage:
        with _SinkVerbosityLevel(2):
            LogVerbose("Test")
    Will only write "Test" in script level (2 above)
    """
    PushVerbosityLevel(GetVerbosityLevel() - n)
    try:
        yield
    except:
        raise
    finally:
        PopVerbosityLevel()


def _QualifiedAtomNotation(a):
    """Return a parsable string of the atom in the format:
    ChainName.ResidueNumber.InsertionCode.AtomName."""
    resnum = a.residue.number
    return "{cname}.{rnum}.{ins_code}.{aname}".format(
        cname=a.chain.name,
        rnum=resnum.num,
        ins_code=resnum.ins_code.strip("\u0000"),
        aname=a.name,
    )


def _QualifiedResidueNotation(r):
    """Return a parsable string of the residue in the format:
    ChainName.ResidueNumber.InsertionCode."""
    resnum = r.number
    return "{cname}.{rnum}.{ins_code}".format(
        cname=r.chain.name,
        rnum=resnum.num,
        ins_code=resnum.ins_code.strip("\u0000"),
    )

class LigandScorer:
    """ Scorer to compute various small molecule ligand (non polymer) scores.

    :class:`LigandScorer` is an abstract base class dealing with all the setup,
    data storage, enumerating ligand symmetries and target/model ligand
    matching/assignment. But actual score computation is delegated to child
    classes.

    At the moment, two such classes are available:

    * :class:`ost.mol.alg.ligand_scoring_lddtpli.LDDTPLIScorer`
      that assesses the conservation of protein-ligand
      contacts (LDDT-PLI);
    * :class:`ost.mol.alg.ligand_scoring_scrmsd.SCRMSDScorer`
      that computes a binding-site superposed, symmetry-corrected RMSD
      (BiSyRMSD) and ligand pocket LDDT (LDDT-LP).

    All versus all scores are available through the lazily computed
    :attr:`score_matrix`. However, many things can go wrong... be it even
    something as simple as two ligands not matching. Error states therefore
    encode scoring issues. An Issue for a particular ligand is indicated by a
    non-zero state in :attr:`model_ligand_states`/:attr:`target_ligand_states`.
    This invalidates pairwise scores of such a ligand with all other ligands.
    This and other issues in pairwise score computation are reported in
    :attr:`state_matrix` which has the same size as :attr:`score_matrix`.
    Only if the respective location is 0, a valid pairwise score can be
    expected. The states and their meaning can be explored with code::

      for state_code, (short_desc, desc) in scorer_obj.state_decoding.items():
          print(state_code)
          print(short_desc)
          print(desc)

    A common use case is to derive a one-to-one mapping between ligands in
    the model and the target for which :class:`LigandScorer` provides an
    automated :attr:`assignment` procedure.
    By default, only exact matches between target and model ligands are
    considered. This is a problem when the target only contains a subset
    of the expected atoms (for instance if atoms are missing in an
    experimental structure, which often happens in the PDB). With
    `substructure_match=True`, complete model ligands can be scored against
    partial target ligands. One problem with this approach is that it is
    very easy to find good matches to small, irrelevant ligands like EDO, CO2
    or GOL. The assignment algorithm therefore considers the coverage,
    expressed as the fraction of atoms of the model ligand atoms covered in the
    target. Higher coverage matches are prioritized, but a match with a better
    score will be preferred if it falls within a window of `coverage_delta`
    (by default 0.2) of a worse-scoring match. As a result, for instance,
    with a delta of 0.2, a low-score match with coverage 0.96 would be
    preferred over a high-score match with coverage 0.70.

    Assumptions:

    Unlike most of OpenStructure, this class does not assume that the ligands
    (either for the model or the target) are part of the PDB component
    dictionary. They may have arbitrary residue names. Residue names do not
    have to match between the model and the target. Matching is based on
    the calculation of isomorphisms which depend on the atom element name and
    atom connectivity (bond order is ignored).
    It is up to the caller to ensure that the connectivity of atoms is properly
    set before passing any ligands to this class. Ligands with improper
    connectivity will lead to bogus results.

    This only applies to the ligand. The rest of the model and target
    structures (protein, nucleic acids) must still follow the usual rules and
    contain only residues from the compound library. Structures are cleaned up
    according to constructor documentation. We advise to
    use the :func:`ost.mol.alg.scoring_base.MMCIFPrep` and
    :func:`ost.mol.alg.scoring_base.PDBPrep` for loading which already
    clean hydrogens and, in the case of MMCIF, optionally extract ligands ready
    to be used by the :class:`LigandScorer` based on "non-polymer" entity types.
    In case of PDB file format, ligands must be loaded separately as SDF files.

    Only polymers (protein and nucleic acids) of model and target are considered
    for ligand binding sites. The
    :class:`ost.mol.alg.chain_mapping.ChainMapper` is used to enumerate possible
    mappings of these chains. In short: identical chains in the target are
    grouped based on pairwise sequence identity
    (see pep_seqid_thr/nuc_seqid_thr param). Each model chain is assigned to
    one of these groups (see mdl_map_pep_seqid_thr/mdl_map_nuc_seqid_thr param).
    To avoid spurious matches, only polymers of a certain length are considered
    in this matching procedure (see min_pep_length/min_nuc_length param).
    Shorter polymers are never mapped and do not contribute to scoring.

    Here is an example of how to setup a scorer::

        from ost.mol.alg.ligand_scoring_scrmsd import SCRMSDScorer
        from ost.mol.alg.scoring_base import MMCIFPrep
        from ost.mol.alg.scoring_base import PDBPrep

        # Load data
        # Structure model in PDB format, containing the receptor only
        model = PDBPrep("path_to_model.pdb")
        # Ligand model as SDF file
        model_ligand = io.LoadEntity("path_to_ligand.sdf", format="sdf")
        # Target loaded from mmCIF, containing the ligand
        target, target_ligands = MMCIFPrep("path_to_target.cif",
                                           extract_nonpoly=True)

        # Setup scorer object and compute SCRMSD
        model_ligands = [model_ligand.Select("ele != H")]
        sc = SCRMSDScorer(model, target, model_ligands, target_ligands)

        # Perform assignment and read respective scores
        for lig_pair in sc.assignment:
            trg_lig = sc.target_ligands[lig_pair[0]]
            mdl_lig = sc.model_ligands[lig_pair[1]]
            score = sc.score_matrix[lig_pair[0], lig_pair[1]]
            print(f"Score for {trg_lig} and {mdl_lig}: {score}")

        # check cleanup in model and target structure:
        print("model cleanup:", sc.model_cleanup_log)
        print("target cleanup:", sc.target_cleanup_log)


    :param model: Model structure - a deep copy is available as :attr:`model`.
                  The model undergoes the following cleanup steps which are
                  dependent on :class:`ost.conop.CompoundLib` returned by
                  :func:`ost.conop.GetDefaultLib`: 1) removal
                  of hydrogens, 2) removal of residues for which there is no
                  entry in :class:`ost.conop.CompoundLib`, 3) removal of
                  residues that are not peptide linking or nucleotide linking
                  according to :class:`ost.conop.CompoundLib` 4) removal of
                  atoms that are not defined for respective residues in
                  :class:`ost.conop.CompoundLib`. Except step 1), every cleanup
                  is logged with :class:`ost.LogLevel` Warning and a report is
                  available as :attr:`model_cleanup_log`.
    :type model: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param target: Target structure - same processing as *model*. 
    :type target: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param model_ligands: Model ligands, as a list of
                          :class:`ost.mol.ResidueHandle`/
                          :class:`ost.mol.ResidueView`/
                          :class:`ost.mol.EntityHandle`/
                          :class:`ost.mol.EntityView`. For
                          :class:`ost.mol.EntityHandle`/
                          :class:`ost.mol.EntityView`, each residue is
                          considered to be an individual ligand.
                          All ligands are copied into a separate
                          :class:`ost.mol.EntityHandle` available as
                          :attr:`model_ligand_ent` and the respective
                          list of ligands is available as :attr:`model_ligands`.
    :type model_ligands: :class:`list`
    :param target_ligands: Target ligands, same processing as model ligands. 
    :type target_ligands: :class:`list`
    :param resnum_alignments: Whether alignments between chemically equivalent
                              chains in *model* and *target* can be computed
                              based on residue numbers. This can be assumed in
                              benchmarking setups such as CAMEO/CASP.
    :type resnum_alignments: :class:`bool`
    :param substructure_match: Set this to True to allow incomplete (i.e.
                               partially resolved) target ligands.
    :type substructure_match: :class:`bool`
    :param coverage_delta: the coverage delta for partial ligand assignment.
    :type coverage_delta: :class:`float`
    :param max_symmetries: If more than that many isomorphisms exist for
                           a target-ligand pair, it will be ignored and reported
                           as unassigned.
    :type max_symmetries: :class:`int`
    :param min_pep_length: Relevant parameter if short peptides are involved in
                           the polymer binding site. Minimum peptide length for
                           a chain to be considered in chain mapping.
                           The chain mapping algorithm first performs an all vs.
                           all pairwise sequence alignment to identify \"equal\"
                           chains within the target structure. We go for simple
                           sequence identity there. Short sequences can be
                           problematic as they may produce high sequence identity
                           alignments by pure chance.
    :type min_pep_length: :class:`int`
    :param min_nuc_length: Same for nucleotides
    :type min_nuc_length: :class:`int`
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

    def __init__(self, model, target, model_ligands, target_ligands,
                 resnum_alignments=False, substructure_match=False,
                 coverage_delta=0.2, max_symmetries=1e5,
                 rename_ligand_chain=False, min_pep_length = 6,
                 min_nuc_length = 4, pep_seqid_thr = 95.,
                 nuc_seqid_thr = 95.,
                 mdl_map_pep_seqid_thr = 0.,
                 mdl_map_nuc_seqid_thr = 0.,
                 seqres = None,
                 trg_seqres_mapping = None):

        if isinstance(model, mol.EntityView):
            self._model = mol.CreateEntityFromView(model, False)
        elif isinstance(model, mol.EntityHandle):
            self._model = model.Copy()
        else:
            raise RuntimeError("model must be of type EntityView/EntityHandle")

        if isinstance(target, mol.EntityView):
            self._target = mol.CreateEntityFromView(target, False)
        elif isinstance(target, mol.EntityHandle):
            self._target = target.Copy()
        else:
            raise RuntimeError("target must be of type EntityView/EntityHandle")

        clib = conop.GetDefaultLib()
        if not clib:
            ost.LogError("A compound library is required. "
                         "Please refer to the OpenStructure website: "
                         "https://openstructure.org/docs/conop/compoundlib/.")
            raise RuntimeError("No compound library found")
        self._target, self._target_cleanup_log = \
        self._cleanup_polymer_ent(self._target, clib)
        self._model, self._model_cleanup_log = \
        self._cleanup_polymer_ent(self._model, clib)

        # keep ligands separate from polymer entities
        self._target_ligand_ent = mol.CreateEntity()
        self._model_ligand_ent = mol.CreateEntity()
        target_ligand_ent_ed = self._target_ligand_ent.EditXCS(mol.BUFFERED_EDIT)
        model_ligand_ent_ed = self._model_ligand_ent.EditXCS(mol.BUFFERED_EDIT)

        self._target_ligands = list()
        for l in target_ligands:
            if isinstance(l, mol.EntityView) or isinstance(l, mol.EntityHandle):
                for r in l.residues:
                    self._target_ligands.append(self._copy_ligand(r, self._target_ligand_ent,
                                                                  target_ligand_ent_ed,
                                                                  rename_ligand_chain))
            elif isinstance(l, mol.ResidueHandle) or isinstance(l, mol.ResidueView):
                self._target_ligands.append(self._copy_ligand(l, self._target_ligand_ent,
                                                              target_ligand_ent_ed,
                                                              rename_ligand_chain))  
            else:
                raise RuntimeError("ligands must be of type EntityView/"
                                   "EntityHandle/ResidueView/ResidueHandle")
            
        self._model_ligands = list()
        for l in model_ligands:
            if isinstance(l, mol.EntityView) or isinstance(l, mol.EntityHandle):
                for r in l.residues:
                    self._model_ligands.append(self._copy_ligand(r, self._model_ligand_ent,
                                                                 model_ligand_ent_ed,
                                                                 rename_ligand_chain))
            elif isinstance(l, mol.ResidueHandle) or isinstance(l, mol.ResidueView):
                self._model_ligands.append(self._copy_ligand(l, self._model_ligand_ent,
                                                             model_ligand_ent_ed,
                                                             rename_ligand_chain))  
            else:
                raise RuntimeError("ligands must be of type EntityView/"
                                   "EntityHandle/ResidueView/ResidueHandle")


        if len(self.model_ligands) == 0:
            LogWarning("No ligands in the model")
            if len(self.target_ligands) == 0:
                raise ValueError("No ligand in the model and in the target")

        self._resnum_alignments = resnum_alignments
        self._substructure_match = substructure_match
        self._coverage_delta = coverage_delta
        self._max_symmetries = max_symmetries
        self._min_pep_length = min_pep_length
        self._min_nuc_length = min_nuc_length
        self._pep_seqid_thr = pep_seqid_thr
        self._nuc_seqid_thr = nuc_seqid_thr
        self._mdl_map_pep_seqid_thr = mdl_map_pep_seqid_thr
        self._mdl_map_nuc_seqid_thr = mdl_map_nuc_seqid_thr
        self._seqres = seqres
        self._trg_seqres_mapping = trg_seqres_mapping

        # lazily computed attributes
        self.__chain_mapper = None
        self.__chem_mapping = None
        self.__chem_group_alns = None
        self.__ref_mdl_alns = None
        self.__mdl_chains_without_chem_mapping = None
        self.__chain_mapping_mdl = None

        # keep track of states
        # simple integers instead of enums - encoding available in
        # self.state_decoding
        self._state_matrix = None
        self._model_ligand_states = None
        self._target_ligand_states = None

        # score matrices
        self._score_matrix = None
        self._coverage_matrix = None
        self._aux_matrix = None

        # assignment and derived data
        self._assignment = None
        self._score_dict = None
        self._aux_dict = None

        # human readable description of states - child class must extend with
        # with child class specific states
        # each state code comes with a tuple of two elements:
        # 1) short description 2) human readable description
        # The actual states are set in _compute_scores in :class:`LigandScorer`
        # or _compute_score of the child class.
        if self.substructure_match:
            iso = "subgraph isomorphism"
        else:
            iso = "full graph isomorphism"

        self.state_decoding = \
        {0: ("OK", "OK"),
         1: ("identity", f"Ligands could not be matched (by {iso})"),
         2: ("symmetries", "Too many symmetries between ligand atoms were "
             "found - increasing max_symmetries might help"),
         3: ("no_iso", "No full isomorphic match could be found - enabling "
             "substructure_match might allow a match"),
         4: ("disconnected", "Ligand graph is disconnected"),
         5: ("stoichiometry", "Ligand was already assigned to another ligand "
             "(different stoichiometry)"),
         6: ("single_ligand_issue", "Cannot compute valid pairwise score as "
             "either model or target ligand have non-zero state."),
         9: ("unknown", "An unknown error occured in LigandScorer")}


    @property
    def model(self):
        """ Model receptor structure

        Processed according to docs in :class:`LigandScorer` constructor
        """
        return self._model

    @property
    def target(self):
        """ Target receptor structure

        Processed according to docs in :class:`LigandScorer` constructor
        """
        return self._target

    @property
    def model_cleanup_log(self):
        """ Reports residues/atoms that were removed in model during cleanup

        Residues and atoms are described as :class:`str` in format
        <chain_name>.<resnum>.<ins_code> (residue) and
        <chain_name>.<resnum>.<ins_code>.<aname> (atom).

        :class:`dict` with keys:
        
        * 'cleaned_residues': another :class:`dict` with keys:

          * 'no_clib': residues that have been removed because no entry could be
            found in :class:`ost.conop.CompoundLib`
          * 'not_linking': residues that have been removed because they're not
            peptide or nucleotide linking according to
            :class:`ost.conop.CompoundLib`

        * 'cleaned_atoms': another :class:`dict` with keys:

          * 'unknown_atoms': atoms that have been removed as they're not part
            of their respective residue according to
            :class:`ost.conop.CompoundLib`
        """
        return self._model_cleanup_log

    @property
    def target_cleanup_log(self):
        """ Same for target
        """
        return self._target_cleanup_log  

    @property
    def model_ligands(self):
        """ Residues representing model ligands

        :class:`list` of :class:`ost.mol.ResidueHandle`
        """
        return self._model_ligands    

    @property
    def target_ligands(self):
        """ Residues representing target ligands

        :class:`list` of :class:`ost.mol.ResidueHandle`
        """
        return self._target_ligands 

    @property
    def resnum_alignments(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._resnum_alignments

    @property
    def min_pep_length(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._min_pep_length
    
    @property
    def min_nuc_length(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._min_nuc_length

    @property
    def pep_seqid_thr(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._pep_seqid_thr
    
    @property
    def nuc_seqid_thr(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._nuc_seqid_thr

    @property
    def mdl_map_pep_seqid_thr(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._mdl_map_pep_seqid_thr
    
    @property
    def mdl_map_nuc_seqid_thr(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._mdl_map_nuc_seqid_thr

    @property
    def seqres(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._seqres

    @property
    def trg_seqres_mapping(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._trg_seqres_mapping

    @property
    def substructure_match(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._substructure_match

    @property
    def coverage_delta(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._coverage_delta

    @property
    def max_symmetries(self):
        """ Given at :class:`LigandScorer` construction
        """
        return self._max_symmetries

    @property
    def state_matrix(self):
        """ Encodes states of ligand pairs

        Ligand pairs can be matched and a valid score can be expected if
        respective location in this matrix is 0.
        Target ligands are in rows, model ligands in columns. States are encoded
        as integers <= 9. Larger numbers encode errors for child classes.
        Use something like ``self.state_decoding[3]`` to get a decscription.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._state_matrix is None:
            self._compute_scores()
        return self._state_matrix

    @property
    def model_ligand_states(self):
        """ Encodes states of model ligands

        Non-zero state in any of the model ligands invalidates the full
        respective column in :attr:`~state_matrix`.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._model_ligand_states is None:
            self._compute_scores()
        return self._model_ligand_states

    @property
    def target_ligand_states(self):
        """ Encodes states of target ligands

        Non-zero state in any of the target ligands invalidates the full
        respective row in :attr:`~state_matrix`.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._target_ligand_states is None:
            self._compute_scores()
        return self._target_ligand_states

    @property
    def score_matrix(self):
        """ Get the matrix of scores.

        Target ligands are in rows, model ligands in columns.

        NaN values indicate that no value could be computed (i.e. different
        ligands). In other words: values are only valid if the respective
        location in :attr:`~state_matrix` is 0. 

        :rtype: :class:`~numpy.ndarray`
        """
        if self._score_matrix is None:
            self._compute_scores() 
        return self._score_matrix

    @property
    def coverage_matrix(self):
        """ Get the matrix of model ligand atom coverage in the target.

        Target ligands are in rows, model ligands in columns.

        NaN values indicate that no value could be computed (i.e. different
        ligands). In other words: values are only valid if the respective
        location in :attr:`~state_matrix` is 0. If `substructure_match=False`,
        only full match isomorphisms are considered, and therefore only values
        of 1.0 can be observed.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._coverage_matrix is None:
            self._compute_scores()
        return self._coverage_matrix

    @property
    def aux_matrix(self):
        """ Get the matrix of scorer specific auxiliary data.

        Target ligands are in rows, model ligands in columns.

        Auxiliary data consists of arbitrary data dicts which allow a child
        class to provide additional information for a scored ligand pair.
        empty dictionaries indicate that the child class simply didn't return
        anything or that no value could be computed (e.g. different ligands).
        In other words: values are only valid if respective location in the
        :attr:`~state_matrix` is 0.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._aux_matrix is None:
            self._compute_scores()
        return self._aux_matrix

    @property
    def assignment(self):
        """ Ligand assignment based on computed scores

        Implements a greedy algorithm to assign target and model ligands
        with each other. Starts from each valid ligand pair as indicated
        by a state of 0 in :attr:`state_matrix`. Each iteration first selects
        high coverage pairs. Given max_coverage defined as the highest
        coverage observed in the available pairs, all pairs with coverage
        in [max_coverage-*coverage_delta*, max_coverage] are selected.
        The best scoring pair among those is added to the assignment
        and the whole process is repeated until there are no ligands to
        assign anymore.

        :rtype: :class:`list` of :class:`tuple` (trg_lig_idx, mdl_lig_idx)
        """
        if self._assignment is None:
            self._assignment = list()
            # Build working array that contains tuples for all mdl/trg ligand
            # pairs with valid score as indicated by a state of 0:
            # (score, coverage, trg_ligand_idx, mdl_ligand_idx)
            tmp = list()
            for trg_idx in range(self.score_matrix.shape[0]):
                for mdl_idx in range(self.score_matrix.shape[1]):
                    if self.state_matrix[trg_idx, mdl_idx] == 0:
                        tmp.append((self.score_matrix[trg_idx, mdl_idx],
                                    self.coverage_matrix[trg_idx, mdl_idx],
                                    trg_idx, mdl_idx))

            # sort by score, such that best scoring item is in front
            if self._score_dir() == '+':
                tmp.sort(reverse=True)
            elif self._score_dir() == '-':
                tmp.sort()
            else:
                raise RuntimeError("LigandScorer._score_dir must return one in "
                                   "['+', '-']")

            LogScript("Computing ligand assignment")
            while len(tmp) > 0:
                # select high coverage ligand pairs in working array
                coverage_thresh = max([x[1] for x in tmp]) - self.coverage_delta
                top_coverage = [x for x in tmp if x[1] >= coverage_thresh]

                # working array is sorted by score => just pick first one
                a = top_coverage[0][2] # selected trg_ligand_idx
                b = top_coverage[0][3] # selected mdl_ligand_idx
                self._assignment.append((a, b))

                # kick out remaining pairs involving these ligands
                tmp = [x for x in tmp if (x[2] != a and x[3] != b)]

        return self._assignment

    @property
    def score(self):
        """ Get a dictionary of score values, keyed by model ligand

        Extract score with something like:
        ``scorer.score[lig.GetChain().GetName()][lig.GetNumber()]``.
        The returned scores are based on :attr:`~assignment`.

        :rtype: :class:`dict`
        """
        if self._score_dict is None:
            self._score_dict = dict()
            for (trg_lig_idx, mdl_lig_idx) in self.assignment:
                mdl_lig = self.model_ligands[mdl_lig_idx]
                cname = mdl_lig.GetChain().GetName()
                rnum = mdl_lig.GetNumber()
                if cname not in self._score_dict:
                    self._score_dict[cname] = dict()
                score = self.score_matrix[trg_lig_idx, mdl_lig_idx]
                self._score_dict[cname][rnum] = score
        return self._score_dict

    @property
    def aux(self):
        """ Get a dictionary of score details, keyed by model ligand
 
        Extract dict with something like:
        ``scorer.score[lig.GetChain().GetName()][lig.GetNumber()]``.
        The returned info dicts are based on :attr:`~assignment`. The content is
        documented in the respective child class.

        :rtype: :class:`dict`
        """
        if self._aux_dict is None:
            self._aux_dict = dict()
            for (trg_lig_idx, mdl_lig_idx) in self.assignment:
                mdl_lig = self.model_ligands[mdl_lig_idx]
                cname = mdl_lig.GetChain().GetName()
                rnum = mdl_lig.GetNumber()
                if cname not in self._aux_dict:
                    self._aux_dict[cname] = dict()
                d = self.aux_matrix[trg_lig_idx, mdl_lig_idx]
                self._aux_dict[cname][rnum] = d
        return self._aux_dict


    @property
    def unassigned_target_ligands(self):
        """ Get indices of target ligands which are not assigned 

        :rtype: :class:`list` of :class:`int`
        """
        # compute on-the-fly, no need for caching
        assigned = set([x[0] for x in self.assignment])
        return [x for x in range(len(self.target_ligands)) if x not in assigned]

    @property
    def unassigned_model_ligands(self):
        """ Get indices of model ligands which are not assigned 

        :rtype: :class:`list` of :class:`int`
        """
        # compute on-the-fly, no need for caching
        assigned = set([x[1] for x in self.assignment])
        return [x for x in range(len(self.model_ligands)) if x not in assigned]

    def get_target_ligand_state_report(self, trg_lig_idx):
        """ Get summary of states observed with respect to all model ligands

        Mainly for debug purposes 

        :param trg_lig_idx: Index of target ligand for which report should be
                            generated
        :type trg_lig_idx: :class:`int`
        """
        return self._get_report(self.target_ligand_states[trg_lig_idx],
                                self.state_matrix[trg_lig_idx,:])

    def get_model_ligand_state_report(self, mdl_lig_idx):
        """ Get summary of states observed with respect to all target ligands

        Mainly for debug purposes 

        :param mdl_lig_idx: Index of model ligand for which report should be
                            generated
        :type mdl_lig_idx: :class:`int`
        """
        return self._get_report(self.model_ligand_states[mdl_lig_idx],
                                self.state_matrix[:, mdl_lig_idx])

    def _get_report(self, ligand_state, pair_states):
        """ Helper
        """
        pair_report = list()
        for s in np.unique(pair_states):
            desc = self.state_decoding[s]
            indices = np.flatnonzero(pair_states == s).tolist()
            pair_report.append({"state": s,
                                "short desc": desc[0],
                                "desc": desc[1],
                                "indices": indices})

        desc = self.state_decoding[ligand_state]
        ligand_report = {"state": ligand_state,
                         "short desc": desc[0],
                         "desc": desc[1]}

        return (ligand_report, pair_report)

    def guess_target_ligand_unassigned_reason(self, trg_lig_idx):
        """ Makes an educated guess why target ligand is not assigned

        This either returns actual error states or custom states that are
        derived from them. Currently, the following reasons are reported:

        * `no_ligand`: there was no ligand in the model.
        * `disconnected`: the ligand graph was disconnected.
        * `identity`: the ligand was not found in the model (by graph
          isomorphism). Check your ligand connectivity.
        * `no_iso`: no full isomorphic match could be found. Try enabling
          `substructure_match=True` if the target ligand is incomplete.
        * `symmetries`: too many symmetries were found (by graph isomorphisms).
          Try to increase `max_symmetries`.
        * `stoichiometry`: there was a possible assignment in the model, but
          the model ligand was already assigned to a different target ligand.
          This indicates different stoichiometries.
        * `no_contact` (LDDT-PLI only): There were no LDDT contacts between
          the binding site and the ligand, and LDDT-PLI is undefined.
        * `target_binding_site` (SCRMSD only): no polymer residues were in
          proximity of the target ligand.
        * `model_binding_site` (SCRMSD only): the binding site was not found
          in the model. Either the binding site was not modeled or the model
          ligand was positioned too far in combination with
          `full_bs_search=False`.

        :param trg_lig_idx: Index of target ligand
        :type trg_lig_idx: :class:`int`
        :returns: :class:`tuple` with two elements: 1) keyword 2) human readable
                  sentence describing the issue, (\"unknown\",\"unknown\") if
                  nothing obvious can be found.
        :raises: :class:`RuntimeError` if specified target ligand is assigned
        """
        if trg_lig_idx not in self.unassigned_target_ligands:
            raise RuntimeError("Specified target ligand is not unassigned")

        # hardcoded tuple if there is simply nothing we can assign it to
        if len(self.model_ligands) == 0:
            return ("no_ligand", "No ligand in the model")

        # if something with the ligand itself is wrong, we can be pretty sure
        # thats why the ligand is unassigned
        if self.target_ligand_states[trg_lig_idx] != 0:
            return self.state_decoding[self.target_ligand_states[trg_lig_idx]]

        # The next best guess comes from looking at pair states
        tmp = np.unique(self.state_matrix[trg_lig_idx,:])

        # In case of any 0, it could have been assigned so it's probably 
        # just not selected due to different stoichiometry - this is no
        # defined state, we just return a hardcoded tuple in this case
        if 0 in tmp:
            return ("stoichiometry",
                    "Ligand was already assigned to an other "
                    "model ligand (different stoichiometry)")

        # maybe it's a symmetry issue?
        if 2 in tmp:
            return self.state_decoding[2]

        # if the state is 6 (single_ligand_issue), there is an issue with its
        # target counterpart.
        if 6 in tmp:
            mdl_idx = np.where(self.state_matrix[trg_lig_idx,:]==6)[0]
            for i in mdl_idx:
                if self.model_ligand_states[i] == 0:
                    raise RuntimeError("This should never happen")
                if self.model_ligand_states[i] != 4 or len(tmp) == 1:
                    # Don't report disconnected if only 1 model ligand is
                    # disconnected, unless that's the only reason
                    return self.state_decoding[self.model_ligand_states[i]]

        # get rid of remaining single ligand issues (only disconnected error)
        if 6 in tmp and len(tmp) > 1:
            tmp = tmp[tmp!=6]

        # prefer everything over identity state
        if 1 in tmp and len(tmp) > 1:
            tmp = tmp[tmp!=1]

        # just return whatever is left
        return self.state_decoding[tmp[0]]


    def guess_model_ligand_unassigned_reason(self, mdl_lig_idx):
        """ Makes an educated guess why model ligand is not assigned

        This either returns actual error states or custom states that are
        derived from them. Currently, the following reasons are reported:

        * `no_ligand`: there was no ligand in the target.
        * `disconnected`: the ligand graph is disconnected.
        * `identity`: the ligand was not found in the target (by graph or
          subgraph isomorphism). Check your ligand connectivity.
        * `no_iso`: no full isomorphic match could be found. Try enabling
          `substructure_match=True` if the target ligand is incomplete.
        * `symmetries`: too many symmetries were found (by graph isomorphisms).
          Try to increase `max_symmetries`.
        * `stoichiometry`: there was a possible assignment in the target, but
          the model target was already assigned to a different model ligand.
          This indicates different stoichiometries.
        * `no_contact` (LDDT-PLI only): There were no LDDT contacts between
          the binding site and the ligand, and LDDT-PLI is undefined.
        * `target_binding_site` (SCRMSD only): a potential assignment was found
          in the target, but there were no polymer residues in proximity of the
          ligand in the target.
        * `model_binding_site` (SCRMSD only): a potential assignment was
          found in the target, but no binding site was found in the model.
          Either the binding site was not modeled or the model ligand was
          positioned too far in combination with `full_bs_search=False`.

        :param mdl_lig_idx: Index of model ligand
        :type mdl_lig_idx: :class:`int`
        :returns: :class:`tuple` with two elements: 1) keyword 2) human readable
                  sentence describing the issue, (\"unknown\",\"unknown\") if
                  nothing obvious can be found.
        :raises: :class:`RuntimeError` if specified model ligand is assigned
        """
        if mdl_lig_idx not in self.unassigned_model_ligands:
            raise RuntimeError("Specified model ligand is not unassigned")

        # hardcoded tuple if there is simply nothing we can assign it to
        if len(self.target_ligands) == 0:
            return ("no_ligand", "No ligand in the target")

        # if something with the ligand itself is wrong, we can be pretty sure
        # thats why the ligand is unassigned
        if self.model_ligand_states[mdl_lig_idx] != 0:
            return self.state_decoding[self.model_ligand_states[mdl_lig_idx]]

        # The next best guess comes from looking at pair states
        tmp = np.unique(self.state_matrix[:,mdl_lig_idx])

        # In case of any 0, it could have been assigned so it's probably 
        # just not selected due to different stoichiometry - this is no
        # defined state, we just return a hardcoded tuple in this case
        if 0 in tmp:
            return ("stoichiometry",
                    "Ligand was already assigned to an other "
                    "model ligand (different stoichiometry)")

        # maybe its a symmetry issue?
        if 2 in tmp:
            return self.state_decoding[2]

        # if the state is 6 (single_ligand_issue), there is an issue with its
        # target counterpart.
        if 6 in tmp:
            trg_idx = np.where(self.state_matrix[:,mdl_lig_idx]==6)[0]
            for i in trg_idx:
                if self.target_ligand_states[i] == 0:
                    raise RuntimeError("This should never happen")
                if self.target_ligand_states[i] != 4 or len(tmp) == 1:
                    # Don't report disconnected if only 1 target ligand is
                    # disconnected, unless that's the only reason
                    return self.state_decoding[self.target_ligand_states[i]]

        # get rid of remaining single ligand issues (only disconnected error)
        if 6 in tmp and len(tmp) > 1:
            tmp = tmp[tmp!=6]

        # prefer everything over identity state
        if 1 in tmp and len(tmp) > 1:
            tmp = tmp[tmp!=1]

        # just return whatever is left
        return self.state_decoding[tmp[0]]

    @property
    def unassigned_model_ligands_reasons(self):
        return_dict = dict()
        for i in self.unassigned_model_ligands:
            lig = self.model_ligands[i]
            cname = lig.GetChain().GetName()
            rnum = lig.GetNumber()
            if cname not in return_dict:
                return_dict[cname] = dict()
            return_dict[cname][rnum] = \
            self.guess_model_ligand_unassigned_reason(i)[0]
        return return_dict
    
    @property
    def unassigned_target_ligands_reasons(self):
        return_dict = dict()
        for i in self.unassigned_target_ligands:
            lig = self.target_ligands[i]
            cname = lig.GetChain().GetName()
            rnum = lig.GetNumber()
            if cname not in return_dict:
                return_dict[cname] = dict()
            return_dict[cname][rnum] = \
            self.guess_target_ligand_unassigned_reason(i)[0]
        return return_dict

    @property
    def _chain_mapper(self):
        """ Chain mapper object for the given :attr:`target`.

        Can be used by child classes if needed, constructed with
        *resnum_alignments* flag

        :type: :class:`ost.mol.alg.chain_mapping.ChainMapper`
        """
        if self.__chain_mapper is None:
            with _SinkVerbosityLevel():
                self.__chain_mapper = \
                chain_mapping.ChainMapper(self.target,
                                          n_max_naive=1e9,
                                          resnum_alignments=self.resnum_alignments,
                                          min_pep_length=self.min_pep_length,
                                          min_nuc_length=self.min_nuc_length,
                                          pep_seqid_thr=self.pep_seqid_thr,
                                          nuc_seqid_thr=self.nuc_seqid_thr,
                                          mdl_map_pep_seqid_thr=self.mdl_map_pep_seqid_thr,
                                          mdl_map_nuc_seqid_thr=self.mdl_map_nuc_seqid_thr,
                                          seqres = self.seqres,
                                          trg_seqres_mapping = self.trg_seqres_mapping)
        return self.__chain_mapper

    @property
    def _chem_mapping(self):
        if self.__chem_mapping is None:
            self.__chem_mapping, self.__chem_group_alns, \
            self.__mdl_chains_without_chem_mapping, self.__chain_mapping_mdl = \
            self._chain_mapper.GetChemMapping(self.model)
        return self.__chem_mapping

    @property
    def _chem_group_alns(self):
        if self.__chem_group_alns is None:   
            self.__chem_mapping, self.__chem_group_alns, \
            self.__mdl_chains_without_chem_mapping, self.__chain_mapping_mdl = \
            self._chain_mapper.GetChemMapping(self.model)
        return self.__chem_group_alns

    @property
    def _ref_mdl_alns(self):
        if self.__ref_mdl_alns is None:
            self.__ref_mdl_alns = \
            chain_mapping._GetRefMdlAlns(self._chain_mapper.chem_groups,
                                    self._chain_mapper.chem_group_alignments,
                                    self._chem_mapping,
                                    self._chem_group_alns)
        return self.__ref_mdl_alns
  
    @property
    def _chain_mapping_mdl(self):
        if self.__chain_mapping_mdl is None:
            with _SinkVerbosityLevel():
                self.__chem_mapping, self.__chem_group_alns, \
                self.__mdl_chains_without_chem_mapping, self.__chain_mapping_mdl = \
                self._chain_mapper.GetChemMapping(self.model)
        return self.__chain_mapping_mdl

    @property
    def _mdl_chains_without_chem_mapping(self):
        if self.__mdl_chains_without_chem_mapping is None:
            self.__chem_mapping, self.__chem_group_alns, \
            self.__mdl_chains_without_chem_mapping, self.__chain_mapping_mdl = \
            self._chain_mapper.GetChemMapping(self.model)
        return self.__mdl_chains_without_chem_mapping

    def _compute_scores(self):
        """
        Compute score for every possible target-model ligand pair and store the
        result in internal matrices.
        """
        ##############################
        # Create the result matrices #
        ##############################
        shape = (len(self.target_ligands), len(self.model_ligands))
        self._score_matrix = np.full(shape, np.nan, dtype=np.float32)
        self._coverage_matrix = np.full(shape, np.nan, dtype=np.float32)
        self._state_matrix = np.full(shape, 0, dtype=np.int32)
        self._model_ligand_states = np.zeros(len(self.model_ligands))
        self._target_ligand_states = np.zeros(len(self.target_ligands))
        self._aux_matrix = np.empty(shape, dtype=dict)

        # precompute ligand graphs
        target_graphs = \
        [_ResidueToGraph(l, by_atom_index=True) for l in self.target_ligands]
        model_graphs = \
        [_ResidueToGraph(l, by_atom_index=True) for l in self.model_ligands]

        for g_idx, g in enumerate(target_graphs):
            if not networkx.is_connected(g):
                self._target_ligand_states[g_idx] = 4
                self._state_matrix[g_idx,:] = 6
                msg = "Disconnected graph observed for target ligand "
                msg += str(self.target_ligands[g_idx])
                LogWarning(msg)

        for g_idx, g in enumerate(model_graphs):
            if not networkx.is_connected(g):
                self._model_ligand_states[g_idx] = 4
                self._state_matrix[:,g_idx] = 6
                msg = "Disconnected graph observed for model ligand "
                msg += str(self.model_ligands[g_idx])
                LogWarning(msg)


        LogScript("Computing pairwise scores for all %s x %s ligands" % shape)
        for target_id, target_ligand in enumerate(self.target_ligands):
            LogInfo("Analyzing target ligand %s" % target_ligand)

            if self._target_ligand_states[target_id] == 4:
                # Disconnected graph - already updated states and reported
                # to LogVerbose
                continue 

            for model_id, model_ligand in enumerate(self.model_ligands):
                LogInfo("Comparing to model ligand %s" % model_ligand)

                #########################################################
                # Compute symmetries for given target/model ligand pair #
                #########################################################

                if self._model_ligand_states[model_id] == 4:
                    # Disconnected graph - already updated states and reported
                    # to LogVerbose
                    continue 

                try:
                    symmetries = ComputeSymmetries(
                        model_ligand, target_ligand,
                        substructure_match=self.substructure_match,
                        by_atom_index=True,
                        max_symmetries=self.max_symmetries,
                        model_graph = model_graphs[model_id],
                        target_graph = target_graphs[target_id])
                    LogInfo("Ligands %s and %s symmetry match" % (
                        str(model_ligand), str(target_ligand)))
                except NoSymmetryError:
                    # Ligands are different - skip
                    LogInfo("No symmetry between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    self._state_matrix[target_id, model_id] = 1
                    continue
                except TooManySymmetriesError:
                    # Ligands are too symmetrical - skip
                    LogWarning("Too many symmetries between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    self._state_matrix[target_id, model_id] = 2
                    continue
                except NoIsomorphicSymmetryError:
                    # Ligands are different - skip
                    LogInfo("No isomorphic symmetry between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    self._state_matrix[target_id, model_id] = 3
                    continue
                except DisconnectedGraphError:
                    # this should never happen as we guard against
                    # DisconnectedGraphError when precomputing the graph
                    LogError("Disconnected graph observed for %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    # just set both ligand states to 4
                    self._model_ligand_states[model_id] = 4
                    self._model_ligand_states[target_id] = 4
                    self._state_matrix[target_id, model_id] = 6
                    continue

                #####################################################
                # Compute score by calling the child class _compute #
                #####################################################
                score, pair_state, target_ligand_state, model_ligand_state, aux\
                 = self._compute(symmetries, target_ligand, model_ligand)

                ############
                # Finalize #
                ############

                # Ensure that returned states are associated with a
                # description. This is a requirement when subclassing
                # LigandScorer => state_decoding dict from base class must
                # be modified in subclass constructor
                if pair_state not in self.state_decoding or \
                   target_ligand_state not in self.state_decoding or \
                   model_ligand_state not in self.state_decoding:
                    raise RuntimeError(f"Subclass returned state "
                                       f"\"{state}\" for which no "
                                       f"description is available. Point "
                                       f"the developer of the used scorer "
                                       f"to this error message.")

                # if any of the ligand states is non-zero, this must be marked
                # by the scorer subclass by specifying a specific pair state
                if target_ligand_state != 0 and pair_state != 6:
                    raise RuntimeError("Observed non-zero target-ligand "
                                       "state which must trigger certain "
                                       "pair state. Point the developer of "
                                       "the used scorer to this error message")

                if model_ligand_state != 0 and pair_state != 6:
                    raise RuntimeError("Observed non-zero model-ligand "
                                       "state which must trigger certain "
                                       "pair state. Point the developer of "
                                       "the used scorer to this error message")

                self._state_matrix[target_id, model_id] = pair_state
                self._target_ligand_states[target_id] = target_ligand_state
                self._model_ligand_states[model_id] = model_ligand_state
                if pair_state == 0:
                    if score is None or np.isnan(score):
                        raise RuntimeError("LigandScorer returned invalid "
                                           "score despite 0 error state")
                    # it's a valid score!
                    self._score_matrix[target_id, model_id] = score
                    cvg = len(symmetries[0][0]) / len(model_ligand.atoms)
                    self._coverage_matrix[target_id, model_id] = cvg
                    self._aux_matrix[target_id, model_id] = aux

    def _compute(self, symmetries, target_ligand, model_ligand):
        """ Compute score for specified ligand pair - defined by child class

        Raises :class:`NotImplementedError` if not implemented by child class.

        :param symmetries: Defines symmetries between *target_ligand* and
                           *model_ligand*. Return value of
                           :func:`ComputeSymmetries`
        :type symmetries: :class:`list` of :class:`tuple` with two elements
                          each: 1) :class:`list` of atom indices in
                          *target_ligand* 2) :class:`list` of respective atom
                          indices in *model_ligand*
        :param target_ligand: The target ligand
        :type target_ligand: :class:`ost.mol.ResidueHandle` or
                             :class:`ost.mol.ResidueView`
        :param model_ligand: The model ligand
        :type model_ligand: :class:`ost.mol.ResidueHandle` or
                            :class:`ost.mol.ResidueView`

        :returns: A :class:`tuple` with three elements: 1) a score
                  (:class:`float`) 2) state (:class:`int`).
                  3) auxiliary data for this ligand pair (:class:`dict`).
                  If state is 0, the score and auxiliary data will be
                  added to :attr:`~score_matrix` and :attr:`~aux_matrix` as well
                  as the respective value in :attr:`~coverage_matrix`.
                  Returned score must be valid in this case (not None/NaN).
                  Child specific non-zero states must be >= 10.
        """
        raise NotImplementedError("_compute must be implemented by child "
                                  "class")

    def _score_dir(self):
        """ Return direction of score - defined by child class

        Relevant for ligand assignment. Must return a string in ['+', '-'].
        '+' for ascending scores, i.e. higher is better (lddt etc.)
        '-' for descending scores, i.e. lower is better (rmsd etc.)
        """
        raise NotImplementedError("_score_dir must be implemented by child "
                                  "class")

    def _copy_ligand(self, l, ent, ed, rename_ligand_chain):
        """ Copies ligand into entity and returns residue handle
        """
        new_chain = ent.FindChain(l.chain.name)
        if not new_chain.IsValid():
            new_chain = ed.InsertChain(l.chain.name)
            ed.SetChainType(new_chain, mol.ChainType.CHAINTYPE_NON_POLY)
        else:
            # Does a residue with the same name already exist?
            already_exists = new_chain.FindResidue(l.number).IsValid()
            if already_exists:
                if rename_ligand_chain:
                    chain_ext = 2  # Extend the chain name by this
                    while True:
                        new_chain_name = \
                        l.chain.name + "_" + str(chain_ext)
                        new_chain = ent.FindChain(new_chain_name)
                        if new_chain.IsValid():
                            chain_ext += 1
                            continue
                        else:
                            new_chain = \
                            ed.InsertChain(new_chain_name)
                            break
                    LogInfo("Moved ligand residue %s to new chain %s" % (
                            l.qualified_name, new_chain.name))
                else:
                    msg = \
                    "A residue number %s already exists in chain %s" % (
                        l.number, l.chain.name)
                    raise RuntimeError(msg)

        # Add the residue with its original residue number
        new_res = ed.AppendResidue(new_chain, l.name, l.number)
        # Add atoms
        for old_atom in l.atoms:
            ed.InsertAtom(new_res, old_atom.name, old_atom.pos, 
                element=old_atom.element, occupancy=old_atom.occupancy,
                b_factor=old_atom.b_factor, is_hetatm=old_atom.is_hetatom)
        # Add bonds
        for old_atom in l.atoms:
            for old_bond in old_atom.bonds:
                new_first = new_res.FindAtom(old_bond.first.name)
                new_second = new_res.FindAtom(old_bond.second.name)
                ed.Connect(new_first, new_second, old_bond.bond_order)
        new_res.SetIsLigand(True)
        return new_res

    def _cleanup_polymer_ent(self, ent, clib):
        """ In principle molck light but logs LigandScorer specific warnings

        Only to be applied to polymer entity

        1) removes atoms with elements set to H or D (not logged as there is no
           effect on scoring)
        2) removes residues with no entry in component dictionary
        3) removes all residues that are not peptide_liking or
           nucleotide_linking according component dictionary
        4) removes unknown atoms according to component dictionary
        5) reruns processor
        """

        cleanup_log = {"cleaned_residues": {"no_clib": [],
                                            "not_linking": []},
                       "cleaned_atoms": {"unknown_atoms": []}}

        # 1)
        hydrogen_sel = ent.Select("ele == H or ele == D")
        if hydrogen_sel.GetAtomCount() > 0:
            # just do, no logging as it has no effect on scoring
            ent = ost.mol.CreateEntityFromView(ent.Select(
                      "ele != H and ele != D"), include_exlusive_atoms=False)

        # extract residues/atoms for 2), 3) and 4) 
        not_in_clib = list()
        not_linking = list()
        unknown_atom = list()

        for r in ent.residues:
            comp = clib.FindCompound(r.GetName())
            if comp is None:
                not_in_clib.append(r)
                continue
            if not (comp.IsPeptideLinking() or comp.IsNucleotideLinking()):
                not_linking.append(r)
                continue
            comp_anames = set([a.name for a in comp.atom_specs])
            for a in r.atoms:
                if a.name not in comp_anames:
                    unknown_atom.append(a)

        # 2)
        if len(not_in_clib) > 0:
            cleanup_log["cleaned_residues"]["no_clib"] = \
            [_QualifiedResidueNotation(r) for r in not_in_clib]
            msg = ("Expected all residues in receptor structure to be defined in "
                   "compound library but this is not the case. "
                   "Please refer to the OpenStructure website if an updated "
                   "library is sufficient to solve the problem: "
                   "https://openstructure.org/docs/conop/compoundlib/ "
                   "These residues will not be considered for scoring (which "
                   "may also affect pre-processing steps such as alignments): ")
            msg += ','.join(cleanup_log["cleaned_residues"]["no_clib"])
            ost.LogWarning(msg)
            for r in not_in_clib:
                r.SetIntProp("clib", 1)
            ent = ost.mol.CreateEntityFromView(ent.Select("grclib:0!=1"),
                      include_exlusive_atoms=False)

        # 3)
        if len(not_linking) > 0:
            cleanup_log["cleaned_residues"]["not_linking"] = \
            [_QualifiedResidueNotation(r) for r in not_linking]
            msg = ("Expected all residues in receptor structure to be peptide "
                   "linking or nucleotide linking according to the compound "
                   "library but this is not the case. "
                   "Please refer to the OpenStructure website if an updated "
                   "library is sufficient to solve the problem: "
                   "https://openstructure.org/docs/conop/compoundlib/ "
                   "These residues will not be considered for scoring (which "
                   "may also affect pre-processing steps such as alignments): ")
            msg += ','.join(cleanup_log["cleaned_residues"]["not_linking"])
            ost.LogWarning(msg)
            for r in not_linking:
                r.SetIntProp("linking", 1)
            ent = ost.mol.CreateEntityFromView(ent.Select("grlinking:0!=1"),
                      include_exlusive_atoms=False)

        # 4)
        if len(unknown_atom) > 0:
            cleanup_log["cleaned_atoms"]["unknown_atoms"] = \
            [_QualifiedAtomNotation(a) for a in unknown_atom]
            msg = ("Expected atom names according to the compound library but "
                   "this is not the case."
                   "Please refer to the OpenStructure website if an updated "
                   "library is sufficient to solve the problem: "
                   "https://openstructure.org/docs/conop/compoundlib/ "
                   "These atoms will not be considered for scoring: ")
            msg += ','.join(cleanup_log["cleaned_atoms"]["unknown_atoms"])
            ost.LogWarning(msg)
            for a in unknown_atom:
                a.SetIntProp("unknown", 1)
            ent = ost.mol.CreateEntityFromView(ent.Select("gaunknown:0!=1"),
                      include_exlusive_atoms=False)

        # 5)
        processor = conop.RuleBasedProcessor(clib)
        processor.Process(ent)

        return (ent, cleanup_log)


def _ResidueToGraph(residue, by_atom_index=False):
    """Return a NetworkX graph representation of the residue.

    :param residue: the residue from which to derive the graph
    :type residue: :class:`ost.mol.ResidueHandle` or
                   :class:`ost.mol.ResidueView`
    :param by_atom_index: Set this parameter to True if you need the nodes to
                          be labeled by atom index (within the residue).
                          Otherwise, if False, the nodes will be labeled by
                          atom names.
    :type by_atom_index: :class:`bool`
    :rtype: :class:`~networkx.classes.graph.Graph`

    Nodes are labeled with the Atom's uppercase
    :attr:`~ost.mol.AtomHandle.element`.
    """
    nxg = networkx.Graph()

    for atom in residue.atoms:
        nxg.add_node(atom.name, element=atom.element.upper())

    # This will list all edges twice - once for every atom of the pair.
    # But as of NetworkX 3.0 adding the same edge twice has no effect,
    # so we're good.
    nxg.add_edges_from([(
        b.first.name,
        b.second.name) for a in residue.atoms for b in a.GetBondList()])

    if by_atom_index:
        nxg = networkx.relabel_nodes(nxg,
                                     {a: b for a, b in zip(
                                         [a.name for a in residue.atoms],
                                         range(len(residue.atoms)))},
                                     True)
    return nxg

def ComputeSymmetries(model_ligand, target_ligand, substructure_match=False,
                      by_atom_index=False, return_symmetries=True,
                      max_symmetries=1e6, model_graph = None,
                      target_graph = None):
    """Return a list of symmetries (isomorphisms) of the model onto the target
    residues.

    :param model_ligand: The model ligand
    :type model_ligand: :class:`ost.mol.ResidueHandle` or
                        :class:`ost.mol.ResidueView`
    :param target_ligand: The target ligand
    :type target_ligand: :class:`ost.mol.ResidueHandle` or
                         :class:`ost.mol.ResidueView`
    :param substructure_match: Set this to True to allow partial ligands
                               in the reference.
    :type substructure_match: :class:`bool`
    :param by_atom_index: Set this parameter to True if you need the symmetries
                          to refer to atom index (within the residue).
                          Otherwise, if False, the symmetries refer to atom
                          names.
    :type by_atom_index: :class:`bool`
    :type return_symmetries: If Truthy, return the mappings, otherwise simply
                             return True if a mapping is found (and raise if
                             no mapping is found). This is useful to quickly
                             find out if a mapping exist without the expensive
                             step to find all the mappings.
    :type return_symmetries: :class:`bool`
    :param max_symmetries: If more than that many isomorphisms exist, raise
      a :class:`TooManySymmetriesError`. This can only be assessed by
      generating at least that many isomorphisms and can take some time.
    :type max_symmetries: :class:`int`
    :raises: :class:`NoSymmetryError` when no symmetry can be found;
             :class:`NoIsomorphicSymmetryError` in case of isomorphic
             subgraph but *substructure_match* is False;
             :class:`TooManySymmetriesError` when more than `max_symmetries`
             isomorphisms are found; :class:`DisconnectedGraphError` if
             graph for *model_ligand*/*target_ligand* is disconnected.
    """

    # Get the Graphs of the ligands
    if model_graph is None:
        model_graph = _ResidueToGraph(model_ligand,
                                      by_atom_index=by_atom_index)

    if not networkx.is_connected(model_graph):
        msg = "Disconnected graph for model ligand %s" % model_ligand
        raise DisconnectedGraphError(msg)


    if target_graph is None:
        target_graph = _ResidueToGraph(target_ligand,
                                       by_atom_index=by_atom_index)

    if not networkx.is_connected(target_graph):
        msg = "Disconnected graph for target ligand %s" % target_ligand
        raise DisconnectedGraphError(msg)

    # Note the argument order (model, target) which differs from spyrmsd.
    # This is because a subgraph of model is isomorphic to target - but not the
    # opposite as we only consider partial ligands in the reference.
    # Make sure to generate the symmetries correctly in the end
    gm = networkx.algorithms.isomorphism.GraphMatcher(
        model_graph, target_graph, node_match=lambda x, y:
        x["element"] == y["element"])
    if gm.is_isomorphic():
        if not return_symmetries:
            return True
        symmetries = []
        for i, isomorphism in enumerate(gm.isomorphisms_iter()):
            if i >= max_symmetries:
                raise TooManySymmetriesError(
                    "Too many symmetries between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
            symmetries.append((list(isomorphism.values()),
                               list(isomorphism.keys())))
        assert len(symmetries) > 0
        LogVerbose("Found %s isomorphic mappings (symmetries)" % len(symmetries))
    elif gm.subgraph_is_isomorphic() and substructure_match:
        if not return_symmetries:
            return True
        symmetries = []
        for i, isomorphism in enumerate(gm.subgraph_isomorphisms_iter()):
            if i >= max_symmetries:
                raise TooManySymmetriesError(
                    "Too many symmetries between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
            symmetries.append((list(isomorphism.values()),
                               list(isomorphism.keys())))
        assert len(symmetries) > 0
        # Assert that all the atoms in the target are part of the substructure
        assert len(symmetries[0][0]) == len(target_ligand.atoms)
        n_sym = len(symmetries)
        LogVerbose("Found %s subgraph isomorphisms (symmetries)" % n_sym)
    elif gm.subgraph_is_isomorphic():
        LogVerbose("Found subgraph isomorphisms (symmetries), but"
                 " ignoring because substructure_match=False")
        raise NoIsomorphicSymmetryError("No symmetry between %s and %s" % (
            str(model_ligand), str(target_ligand)))
    else:
        LogVerbose("Found no isomorphic mappings (symmetries)")
        raise NoSymmetryError("No symmetry between %s and %s" % (
            str(model_ligand), str(target_ligand)))

    return symmetries

class NoSymmetryError(ValueError):
    """ Exception raised when no symmetry can be found.
    """
    pass

class NoIsomorphicSymmetryError(ValueError):
    """ Exception raised when no isomorphic symmetry can be found

    There would be isomorphic subgraphs for which symmetries can
    be found, but substructure_match is disabled
    """
    pass

class TooManySymmetriesError(ValueError):
    """ Exception raised when too many symmetries are found.
    """
    pass

class DisconnectedGraphError(Exception):
    """ Exception raised when the ligand graph is disconnected.
    """
    pass

# specify public interface
__all__ = ('LigandScorer', 'ComputeSymmetries', 'NoSymmetryError',
           'NoIsomorphicSymmetryError', 'TooManySymmetriesError',
           'DisconnectedGraphError')
