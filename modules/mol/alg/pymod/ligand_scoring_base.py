import numpy as np
import networkx

from ost import mol
from ost import LogWarning, LogScript, LogVerbose, LogDebug
from ost.mol.alg import chain_mapping

class LigandScorer:

    def __init__(self, model, target, model_ligands=None, target_ligands=None,
                 resnum_alignments=False, rename_ligand_chain=False,
                 substructure_match=False, coverage_delta=0.2,
                 max_symmetries=1e5):

        if isinstance(model, mol.EntityView):
            self.model = mol.CreateEntityFromView(model, False)
        elif isinstance(model, mol.EntityHandle):
            self.model = model.Copy()
        else:
            raise RuntimeError("model must be of type EntityView/EntityHandle")

        if isinstance(target, mol.EntityView):
            self.target = mol.CreateEntityFromView(target, False)
        elif isinstance(target, mol.EntityHandle):
            self.target = target.Copy()
        else:
            raise RuntimeError("target must be of type EntityView/EntityHandle")

        # Extract ligands from target
        if target_ligands is None:
            self.target_ligands = self._extract_ligands(self.target)
        else:
            self.target_ligands = self._prepare_ligands(self.target, target,
                                                        target_ligands,
                                                        rename_ligand_chain)
        if len(self.target_ligands) == 0:
            LogWarning("No ligands in the target")

        # Extract ligands from model
        if model_ligands is None:
            self.model_ligands = self._extract_ligands(self.model)
        else:
            self.model_ligands = self._prepare_ligands(self.model, model,
                                                       model_ligands,
                                                       rename_ligand_chain)
        if len(self.model_ligands) == 0:
            LogWarning("No ligands in the model")
            if len(self.target_ligands) == 0:
                raise ValueError("No ligand in the model and in the target")

        self.resnum_alignments = resnum_alignments
        self.rename_ligand_chain = rename_ligand_chain
        self.substructure_match = substructure_match
        self.coverage_delta = coverage_delta
        self.max_symmetries = max_symmetries

        # lazily computed attributes
        self._chain_mapper = None

        # keep track of states
        # simple integers instead of enums - documentation of property describes
        # encoding
        self._states = None

        # score matrices
        self._score_matrix = None
        self._coverage_matrix = None
        self._aux_data = None

    @property
    def states(self):
        """ Encodes states of ligand pairs

        Expect a valid score if respective location in this matrix is 0.
        Target ligands are in rows, model ligands in columns. States are encoded
        as integers <= 9. Larger numbers encode errors for child classes.
        
        * -1: Unknown Error - cannot be matched
        * 0: Ligand pair can be matched and valid score is computed.
        * 1: Ligand pair has no valid symmetry - cannot be matched.
        * 2: Ligand pair has too many symmetries - cannot be matched. 
             You might be able to get a match by increasing *max_symmetries*.
        * 3: Ligand pair has no isomorphic symmetries - cannot be matched.
             Target ligand is subgraph of model ligand. This error only occurs
             if *substructure_match* is False. These cases will likely become
             0 if this flag is enabled.
        * 4: Disconnected graph error - cannot be matched.
             Either target ligand or model ligand has disconnected graph.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._states is None:
            self._compute_scores()
        return self._states

    @property
    def score_matrix(self):
        """ Get the matrix of scores.

        Target ligands are in rows, model ligands in columns.

        NaN values indicate that no value could be computed (i.e. different
        ligands). In other words: values are only valid if respective location
        :attr:`~states` is 0. 

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
        ligands). In other words: values are only valid if respective location
        :attr:`~states` is 0. If `substructure_match=False`, only full
        match isomorphisms are considered, and therefore only values of 1.0
        can be observed.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._coverage_matrix is None:
            self._compute_scores()
        return self._coverage_matrix

    @property
    def aux_data(self):
        """ Get the matrix of scorer specific auxiliary data.

        Target ligands are in rows, model ligands in columns.

        Auxiliary data consists of arbitrary data dicts which allow a child
        class to provide additional information for a scored ligand pair.
        empty dictionaries indicate that no value could be computed
        (i.e. different ligands). In other words: values are only valid if
        respective location :attr:`~states` is 0.

        :rtype: :class:`~numpy.ndarray`
        """
        if self._aux_matrix is None:
            self._compute_scores()
        return self._aux_matrix

    @property
    def chain_mapper(self):
        """ Chain mapper object for the given :attr:`target`.

        Can be used by child classes if needed, constructed with
        *resnum_alignments* flag

        :type: :class:`ost.mol.alg.chain_mapping.ChainMapper`
        """
        if self._chain_mapper is None:
            self._chain_mapper = \
            chain_mapping.ChainMapper(self.target,
                                      n_max_naive=1e9,
                                      resnum_alignments=self.resnum_alignments)
        return self._chain_mapper

    @staticmethod
    def _extract_ligands(entity):
        """Extract ligands from entity. Return a list of residues.

        Assumes that ligands have the :attr:`~ost.mol.ResidueHandle.is_ligand`
        flag set. This is typically the case for entities loaded from mmCIF
        (tested with mmCIF files from the PDB and SWISS-MODEL).
        Legacy PDB files must contain `HET` headers (which is usually the
        case for files downloaded from the PDB but not elsewhere).

        This function performs basic checks to ensure that the residues in this
        chain are not forming polymer bonds (ie peptide/nucleotide ligands) and
        will raise a RuntimeError if this assumption is broken.

        :param entity: the entity to extract ligands from
        :type entity: :class:`~ost.mol.EntityHandle`
        :rtype: :class:`list` of :class:`~ost.mol.ResidueHandle`

        """
        extracted_ligands = []
        for residue in entity.residues:
            if residue.is_ligand:
                if mol.InSequence(residue, residue.next):
                    raise RuntimeError("Residue %s connected in polymer sequen"
                                       "ce %s" % (residue.qualified_name))
                extracted_ligands.append(residue)
                LogVerbose("Detected residue %s as ligand" % residue)
        return extracted_ligands

    @staticmethod
    def _prepare_ligands(new_entity, old_entity, ligands, rename_chain):
        """Prepare the ligands given into a list of ResidueHandles which are
        part of the copied entity, suitable for the model_ligands and
        target_ligands properties.

        This function takes a list of ligands as (Entity|Residue)(Handle|View).
        Entities can contain multiple ligands, which will be considered as
        separate ligands.

        Ligands which are part of the entity are simply fetched in the new
        copied entity. Otherwise, they are copied over to the copied entity.
        """
        extracted_ligands = []

        next_chain_num = 1
        new_editor = None

        def _copy_residue(residue, rename_chain):
            """ Copy the residue into the new chain.
            Return the new residue handle."""
            nonlocal next_chain_num, new_editor

            # Instantiate the editor
            if new_editor is None:
                new_editor = new_entity.EditXCS()

            new_chain = new_entity.FindChain(residue.chain.name)
            if not new_chain.IsValid():
                new_chain = new_editor.InsertChain(residue.chain.name)
            else:
                # Does a residue with the same name already exist?
                already_exists = new_chain.FindResidue(residue.number).IsValid()
                if already_exists:
                    if rename_chain:
                        chain_ext = 2  # Extend the chain name by this
                        while True:
                            new_chain_name = residue.chain.name + "_" + str(chain_ext)
                            new_chain = new_entity.FindChain(new_chain_name)
                            if new_chain.IsValid():
                                chain_ext += 1
                                continue
                            else:
                                new_chain = new_editor.InsertChain(new_chain_name)
                                break
                        LogScript("Moved ligand residue %s to new chain %s" % (
                            residue.qualified_name, new_chain.name))
                    else:
                        msg = "A residue number %s already exists in chain %s" % (
                            residue.number, residue.chain.name)
                        raise RuntimeError(msg)

            # Add the residue with its original residue number
            new_res = new_editor.AppendResidue(new_chain, residue.name, residue.number)
            # Add atoms
            for old_atom in residue.atoms:
                new_editor.InsertAtom(new_res, old_atom.name, old_atom.pos, 
                    element=old_atom.element, occupancy=old_atom.occupancy,
                    b_factor=old_atom.b_factor, is_hetatm=old_atom.is_hetatom)
            # Add bonds
            for old_atom in residue.atoms:
                for old_bond in old_atom.bonds:
                    new_first = new_res.FindAtom(old_bond.first.name)
                    new_second = new_res.FindAtom(old_bond.second.name)
                    new_editor.Connect(new_first, new_second)
            return new_res

        def _process_ligand_residue(res, rename_chain):
            """Copy or fetch the residue. Return the residue handle."""
            new_res = None
            if res.entity.handle == old_entity.handle:
                # Residue is part of the old_entity handle.
                # However, it may not be in the copied one, for instance it may have been a view
                # We try to grab it first, otherwise we copy it
                new_res = new_entity.FindResidue(res.chain.name, res.number)
            if new_res and new_res.valid:
                LogVerbose("Ligand residue %s already in entity" % res.handle.qualified_name)
            else:
                # Residue is not part of the entity, need to copy it first
                new_res = _copy_residue(res, rename_chain)
                LogVerbose("Copied ligand residue %s" % res.handle.qualified_name)
            new_res.SetIsLigand(True)
            return new_res

        for ligand in ligands:
            if isinstance(ligand, mol.EntityHandle) or isinstance(ligand, mol.EntityView):
                for residue in ligand.residues:
                    new_residue = _process_ligand_residue(residue, rename_chain)
                    extracted_ligands.append(new_residue)
            elif isinstance(ligand, mol.ResidueHandle) or isinstance(ligand, mol.ResidueView):
                new_residue = _process_ligand_residue(ligand, rename_chain)
                extracted_ligands.append(new_residue)
            else:
                raise RuntimeError("Ligands should be given as Entity or Residue")

        if new_editor is not None:
            new_editor.UpdateICS()
        return extracted_ligands

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
        self._states = np.full(shape, -1, dtype=np.int32)
        self._aux_data = np.empty(shape, dtype=dict)

        for target_id, target_ligand in enumerate(self.target_ligands):
            LogVerbose("Analyzing target ligand %s" % target_ligand)
            for model_id, model_ligand in enumerate(self.model_ligands):
                LogVerbose("Compare to model ligand %s" % model_ligand)

                #########################################################
                # Compute symmetries for given target/model ligand pair #
                #########################################################
                try:
                    symmetries = ComputeSymmetries(
                        model_ligand, target_ligand,
                        substructure_match=self.substructure_match,
                        by_atom_index=True,
                        max_symmetries=self.max_symmetries)
                    LogVerbose("Ligands %s and %s symmetry match" % (
                        str(model_ligand), str(target_ligand)))
                except NoSymmetryError:
                    # Ligands are different - skip
                    LogVerbose("No symmetry between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    self._states[target_id, model_id] = 1
                    continue
                except TooManySymmetriesError:
                    # Ligands are too symmetrical - skip
                    LogVerbose("Too many symmetries between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    self._states[target_id, model_id] = 2
                    continue
                except NoIsomorphicSymmetryError:
                    # Ligands are different - skip
                    LogVerbose("No isomorphic symmetry between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    self._states[target_id, model_id] = 3
                    continue
                except DisconnectedGraphError:
                    LogVerbose("Disconnected graph observed for %s and %s" % (
                        str(model_ligand), str(target_ligand)))
                    self._states[target_id, model_id] = 4
                    continue

                #####################################################
                # Compute score by calling the child class _compute #
                #####################################################
                score, state, aux = self._compute(symmetries, target_ligand,
                                                  model_ligand)

                ############
                # Finalize #
                ############
                if state != 0:
                    # non-zero error states up to 4 are reserved for base class
                    if state <= 9:
                        raise RuntimeError("Child returned reserved err. state")

                self._states[target_id, model_id] = state
                if state == 0:
                    if score is None or np.isnan(score):
                        raise RuntimeError("LigandScorer returned invalid "
                                           "score despite 0 error state")
                    # it's a valid score!
                    self._score_matrix[target_id, model_id] = score
                    cvg = len(symmetries[0][0]) / len(model_ligand.atoms)
                    self._coverage_matrix[target_id, model_id] = cvg
                    self._aux_data[target_id, model_id] = aux

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
                  added to :attr:`~score_matrix` and :attr:`~aux_data` as well
                  as the respective value in :attr:`~coverage_matrix`.
                  Returned score must be valid in this case (not None/NaN).
                  Child specific non-zero states must be >= 10.
        """
        raise NotImplementedError("_compute must be implemented by child class")


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

    Nodes are labeled with the Atom's uppercase :attr:`~ost.mol.AtomHandle.element`.
    """
    nxg = networkx.Graph()

    for atom in residue.atoms:
        nxg.add_node(atom.name, element=atom.element.upper())

    # This will list all edges twice - once for every atom of the pair.
    # But as of NetworkX 3.0 adding the same edge twice has no effect, so we're good.
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
                      max_symmetries=1e6):
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
             subgraph but *substructure_match* is False.
             :class:`TooManySymmetriesError` when more than `max_symmetries`
             isomorphisms are found.
    """

    # Get the Graphs of the ligands
    model_graph = _ResidueToGraph(model_ligand, by_atom_index=by_atom_index)
    target_graph = _ResidueToGraph(target_ligand, by_atom_index=by_atom_index)

    if not networkx.is_connected(model_graph):
        raise DisconnectedGraphError("Disconnected graph for model ligand %s" % model_ligand)
    if not networkx.is_connected(target_graph):
        raise DisconnectedGraphError("Disconnected graph for target ligand %s" % target_ligand)

    # Note the argument order (model, target) which differs from spyrmsd.
    # This is because a subgraph of model is isomorphic to target - but not the opposite
    # as we only consider partial ligands in the reference.
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
            symmetries.append((list(isomorphism.values()), list(isomorphism.keys())))
        assert len(symmetries) > 0
        LogDebug("Found %s isomorphic mappings (symmetries)" % len(symmetries))
    elif gm.subgraph_is_isomorphic() and substructure_match:
        if not return_symmetries:
            return True
        symmetries = []
        for i, isomorphism in enumerate(gm.subgraph_isomorphisms_iter()):
            if i >= max_symmetries:
                raise TooManySymmetriesError(
                    "Too many symmetries between %s and %s" % (
                        str(model_ligand), str(target_ligand)))
            symmetries.append((list(isomorphism.values()), list(isomorphism.keys())))
        assert len(symmetries) > 0
        # Assert that all the atoms in the target are part of the substructure
        assert len(symmetries[0][0]) == len(target_ligand.atoms)
        LogDebug("Found %s subgraph isomorphisms (symmetries)" % len(symmetries))
    elif gm.subgraph_is_isomorphic():
        LogDebug("Found subgraph isomorphisms (symmetries), but"
                 " ignoring because substructure_match=False")
        raise NoIsomorphicSymmetryError("No symmetry between %s and %s" % (
            str(model_ligand), str(target_ligand)))
    else:
        LogDebug("Found no isomorphic mappings (symmetries)")
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
