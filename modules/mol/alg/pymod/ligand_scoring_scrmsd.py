import numpy as np

from ost import LogWarning, LogScript, LogInfo, LogVerbose
from ost import geom
from ost import mol

from ost.mol.alg import ligand_scoring_base


class SCRMSDScorer(ligand_scoring_base.LigandScorer):
    """ :class:`LigandScorer` implementing symmetry corrected RMSD (BiSyRMSD).

    :class:`SCRMSDScorer` computes a score for a specific pair of target/model
    ligands.

    The returned RMSD is based on a binding site superposition.
    The binding site of the target structure is defined as all residues with at
    least one atom within `bs_radius` around the target ligand.
    It only contains protein and nucleic acid residues from chains that
    pass the criteria for the
    :class:`chain mapping <ost.mol.alg.chain_mapping>`. This means ignoring
    other ligands, waters, short polymers as well as any incorrectly connected
    chains that may be in proximity.
    The respective model binding site for superposition is identified by
    naively enumerating all possible mappings of model chains onto their
    chemically equivalent target counterparts from the target binding site.
    The `binding_sites_topn` with respect to lDDT score are evaluated and 
    an RMSD is computed.
    You can either try to map ALL model chains onto the target binding site by
    enabling `full_bs_search` or restrict the model chains for a specific
    target/model ligand pair to the chains with at least one atom within
    *model_bs_radius* around the model ligand. The latter can be significantly
    faster in case of large complexes.
    Symmetry correction is achieved by simply computing an RMSD value for
    each symmetry, i.e. atom-atom assignments of the ligand as given by
    :class:`LigandScorer`. The lowest RMSD value is returned.

    Populates :attr:`LigandScorer.aux_data` with following :class:`dict` keys:

    * rmsd: The BiSyRMSD score
    * lddt_lp: lDDT of the binding pocket used for superposition (lDDT-LP)
    * bs_ref_res: :class:`list` of binding site residues in target
    * bs_ref_res_mapped: :class:`list` of target binding site residues that
      are mapped to model
    * bs_mdl_res_mapped: :class:`list` of same length with respective model
      residues
    * bb_rmsd: Backbone RMSD (CA, C3' for nucleotides; full backbone for 
      binding sites with fewer than 3 residues) for mapped binding site
      residues after superposition
    * target_ligand: The actual target ligand for which the score was computed
    * model_ligand: The actual model ligand for which the score was computed
    * chain_mapping: :class:`dict` with a chain mapping of chains involved in
      binding site - key: trg chain name, value: mdl chain name
    * transform: :class:`geom.Mat4` to transform model binding site onto target
      binding site
    * inconsistent_residues: :class:`list` of :class:`tuple` representing
      residues with inconsistent residue names upon mapping (which is given by
      bs_ref_res_mapped and bs_mdl_res_mapped). Tuples have two elements:
      1) trg residue 2) mdl residue

    :param model: Passed to parent constructor - see :class:`LigandScorer`.
    :type model: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param target: Passed to parent constructor - see :class:`LigandScorer`.
    :type target: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param model_ligands: Passed to parent constructor - see
                          :class:`LigandScorer`.
    :type model_ligands: :class:`list`
    :param target_ligands: Passed to parent constructor - see
                           :class:`LigandScorer`.
    :type target_ligands: :class:`list`
    :param resnum_alignments: Passed to parent constructor - see
                              :class:`LigandScorer`.
    :type resnum_alignments: :class:`bool`
    :param rename_ligand_chain: Passed to parent constructor - see
                                :class:`LigandScorer`.
    :type rename_ligand_chain: :class:`bool`
    :param substructure_match: Passed to parent constructor - see
                               :class:`LigandScorer`.
    :type substructure_match: :class:`bool`
    :param coverage_delta: Passed to parent constructor - see
                           :class:`LigandScorer`.
    :type coverage_delta: :class:`float`
    :param max_symmetries: Passed to parent constructor - see
                           :class:`LigandScorer`.
    :type max_symmetries: :class:`int`
    :param bs_radius: Inclusion radius for the binding site. Residues with
                   atoms within this distance of the ligand will be considered
                   for inclusion in the binding site.
    :type bs_radius: :class:`float`
    :param lddt_lp_radius: lDDT inclusion radius for lDDT-LP.
    :type lddt_lp_radius: :class:`float`
    :param model_bs_radius: inclusion radius for model binding sites.
                            Only used when full_bs_search=False, otherwise the
                            radius is effectively infinite. Only chains with
                            atoms within this distance of a model ligand will
                            be considered in the chain mapping.
    :type model_bs_radius: :class:`float`
    :param binding_sites_topn: maximum number of model binding site
                               representations to assess per target binding
                               site.
    :type binding_sites_topn: :class:`int`
    :param full_bs_search: If True, all potential binding sites in the model
                           are searched for each target binding site. If False,
                           the search space in the model is reduced to chains
                           around (`model_bs_radius` Å) model ligands.
                           This speeds up computations, but may result in
                           ligands not being scored if the predicted ligand
                           pose is too far from the actual binding site.
    :type full_bs_search: :class:`bool`
    :param min_pep_length: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`.
    :type min_pep_length: :class:`int`
    :param min_nuc_length: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`
    :type min_nuc_length: :class:`int`
    :param pep_seqid_thr: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`
    :type pep_seqid_thr: :class:`float`
    :param nuc_seqid_thr: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`
    :type nuc_seqid_thr: :class:`float`
    :param mdl_map_pep_seqid_thr: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`
    :type mdl_map_pep_seqid_thr: :class:`float`
    :param mdl_map_nuc_seqid_thr: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`
    :type mdl_map_nuc_seqid_thr: :class:`float`
    :param seqres: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`
    :type seqres: :class:`ost.seq.SequenceList`
    :param trg_seqres_mapping: See :class:`ost.mol.alg.ligand_scoring_base.LigandScorer`
    :type trg_seqres_mapping: :class:`dict`
    """
    def __init__(self, model, target, model_ligands, target_ligands,
                 resnum_alignments=False, rename_ligand_chain=False,
                 substructure_match=False, coverage_delta=0.2,
                 max_symmetries=1e5, bs_radius=4.0, lddt_lp_radius=15.0,
                 model_bs_radius=25, binding_sites_topn=100000,
                 full_bs_search=False, min_pep_length = 6,
                 min_nuc_length = 4, pep_seqid_thr = 95.,
                 nuc_seqid_thr = 95.,
                 mdl_map_pep_seqid_thr = 0.,
                 mdl_map_nuc_seqid_thr = 0.,
                 seqres=None,
                 trg_seqres_mapping=None):

        super().__init__(model, target, model_ligands, target_ligands,
                         resnum_alignments = resnum_alignments,
                         rename_ligand_chain = rename_ligand_chain,
                         substructure_match = substructure_match,
                         coverage_delta = coverage_delta,
                         max_symmetries = max_symmetries,
                         min_pep_length = min_pep_length,
                         min_nuc_length = min_nuc_length,
                         pep_seqid_thr = pep_seqid_thr,
                         nuc_seqid_thr = nuc_seqid_thr,
                         mdl_map_pep_seqid_thr = mdl_map_pep_seqid_thr,
                         mdl_map_nuc_seqid_thr = mdl_map_nuc_seqid_thr,
                         seqres = seqres,
                         trg_seqres_mapping = trg_seqres_mapping)

        self.bs_radius = bs_radius
        self.lddt_lp_radius = lddt_lp_radius
        self.model_bs_radius = model_bs_radius
        self.binding_sites_topn = binding_sites_topn
        self.full_bs_search = full_bs_search

        # Residues that are in contact with a ligand => binding site
        # defined as all residues with at least one atom within self.radius
        # key: ligand.handle.hash_code, value: EntityView of whatever
        # entity ligand belongs to
        self._binding_sites = dict()

        # cache for GetRepr chain mapping calls
        self._repr = dict()

        # lazily precomputed variables to speedup GetRepr chain mapping calls
        # for localized GetRepr searches
        self._get_repr_input = dict()

        # update state decoding from parent with subclass specific stuff
        self.state_decoding[10] = ("target_binding_site",
                                   "No residues were in proximity of the "
                                   "target ligand.")
        self.state_decoding[11] = ("model_binding_site", "Binding site was not"
                                   " found in the model, i.e. the binding site"
                                   " was not modeled or the model ligand was "
                                   "positioned too far in combination with "
                                   "full_bs_search=False.")
        self.state_decoding[20] = ("unknown",
                                   "Unknown error occured in SCRMSDScorer")

    def _compute(self, symmetries, target_ligand, model_ligand):
        """ Implements interface from parent
        """
        # set default to invalid scores
        best_rmsd_result = {"rmsd": None,
                           "lddt_lp": None,
                           "bs_ref_res": list(),
                           "bs_ref_res_mapped": list(),
                           "bs_mdl_res_mapped": list(),
                           "bb_rmsd": None,
                           "target_ligand": target_ligand,
                           "model_ligand": model_ligand,
                           "chain_mapping": dict(),
                           "transform": geom.Mat4(),
                           "inconsistent_residues": list()}

        representations = self._get_repr(target_ligand, model_ligand)
        # This step can be slow so give some hints in logs
        msg = "Computing BiSyRMSD with %d chain mappings" % len(representations)
        (LogWarning if len(representations) > 10000 else LogInfo)(msg)

        for r in representations:
            rmsd = _SCRMSD_symmetries(symmetries, model_ligand, 
                                      target_ligand, transformation=r.transform)

            if best_rmsd_result["rmsd"] is None or \
               rmsd < best_rmsd_result["rmsd"]:
                best_rmsd_result = {"rmsd": rmsd,
                                    "lddt_lp": r.lDDT,
                                    "bs_ref_res": r.substructure.residues,
                                    "bs_ref_res_mapped": r.ref_residues,
                                    "bs_mdl_res_mapped": r.mdl_residues,
                                    "bb_rmsd": r.bb_rmsd,
                                    "target_ligand": target_ligand,
                                    "model_ligand": model_ligand,
                                    "chain_mapping": r.GetFlatChainMapping(),
                                    "transform": r.transform,
                                    "inconsistent_residues":
                                    r.inconsistent_residues}

        target_ligand_state = 0
        model_ligand_state = 0
        pair_state = 0

        if best_rmsd_result["rmsd"] is not None:
            best_rmsd = best_rmsd_result["rmsd"]
        else:
            # try to identify error states
            best_rmsd = np.nan
            error_state = 20 # unknown error
            N = self._get_target_binding_site(target_ligand).GetResidueCount() 
            if N == 0:
                pair_state = 6 # binding_site
                target_ligand_state = 10
            elif len(representations) == 0:
                pair_state = 11

        return (best_rmsd, pair_state, target_ligand_state, model_ligand_state,
                best_rmsd_result)

    def _score_dir(self):
        """ Implements interface from parent
        """
        return '-'

    def _get_repr(self, target_ligand, model_ligand):

        key = None
        if self.full_bs_search:
            # all possible binding sites, independent from actual model ligand
            key = (target_ligand.handle.hash_code, 0)
        else:
            key = (target_ligand.handle.hash_code,
                   model_ligand.handle.hash_code)

        if key not in self._repr:
            ref_bs = self._get_target_binding_site(target_ligand)
            LogVerbose("%d chains are in proximity of the target ligand: %s" % (
                ref_bs.chain_count, ", ".join([c.name for c in ref_bs.chains])))
            if self.full_bs_search:
                reprs = self._chain_mapper.GetRepr(
                    ref_bs, self.model, inclusion_radius=self.lddt_lp_radius,
                    topn=self.binding_sites_topn)
            else:
                repr_in = self._get_get_repr_input(model_ligand)
                radius = self.lddt_lp_radius
                reprs = self._chain_mapper.GetRepr(ref_bs, self.model,
                                                  inclusion_radius=radius,
                                                  topn=self.binding_sites_topn,
                                                  chem_mapping_result=repr_in)
            self._repr[key] = reprs

        return self._repr[key]

    def _get_target_binding_site(self, target_ligand):

        if target_ligand.handle.hash_code not in self._binding_sites:

            # create view of reference binding site
            ref_residues_hashes = set() # helper to keep track of added residues
            ignored_residue_hashes = {target_ligand.hash_code}
            for ligand_at in target_ligand.atoms:
                close_atoms = self.target.FindWithin(ligand_at.GetPos(),
                                                     self.bs_radius)
                for close_at in close_atoms:
                    # Skip any residue not in the chain mapping target
                    ref_res = close_at.GetResidue()
                    h = ref_res.handle.GetHashCode()
                    if h not in ref_residues_hashes and \
                            h not in ignored_residue_hashes:
                        with ligand_scoring_base._SinkVerbosityLevel(1):
                            view = self._chain_mapper.target.ViewForHandle(ref_res) 
                        if view.IsValid():
                            h = ref_res.handle.GetHashCode()
                            ref_residues_hashes.add(h)
                        elif ref_res.is_ligand:
                            msg = f"Ignoring ligand {ref_res.qualified_name} "
                            msg += "in binding site of "
                            msg += str(target_ligand.qualified_name)
                            LogWarning(msg)
                            ignored_residue_hashes.add(h)
                        elif ref_res.chem_type == mol.ChemType.WATERS:
                            pass # That's ok, no need to warn
                        else:
                            msg = f"Ignoring residue {ref_res.qualified_name} "
                            msg += "in binding site of "
                            msg += str(target_ligand.qualified_name)
                            LogWarning(msg)
                            ignored_residue_hashes.add(h)

            ref_bs = self.target.CreateEmptyView()
            if ref_residues_hashes:
                # reason for doing that separately is to guarantee same ordering
                # of residues as in underlying entity. (Reorder by ResNum seems
                # only available on ChainHandles)
                for ch in self.target.chains:
                    for r in ch.residues:
                        if r.handle.GetHashCode() in ref_residues_hashes:
                            ref_bs.AddResidue(r, mol.ViewAddFlag.INCLUDE_ALL)
                if len(ref_bs.residues) == 0:
                    raise RuntimeError("Failed to add proximity residues to "
                                       "the reference binding site entity")

            self._binding_sites[target_ligand.handle.hash_code] = ref_bs

        return self._binding_sites[target_ligand.handle.hash_code]

    def _get_get_repr_input(self, mdl_ligand):
        if mdl_ligand.handle.hash_code not in self._get_repr_input:

            # figure out what chains in the model are in contact with the ligand
            # that may give a non-zero contribution to lDDT in
            # chain_mapper.GetRepr
            radius = self.model_bs_radius
            chains = set()
            for at in mdl_ligand.atoms:
                with ligand_scoring_base._SinkVerbosityLevel(1):
                    close_atoms = self._chain_mapping_mdl.FindWithin(at.GetPos(),
                                                                     radius)
                for close_at in close_atoms:
                    chains.add(close_at.GetChain().GetName())

            if len(chains) > 0:
                LogVerbose("%d chains are in proximity of the model ligand: %s" % (
                    len(chains), ", ".join(chains)))

                # the chain mapping model which only contains close chains
                query = "cname="
                query += ','.join([mol.QueryQuoteName(x) for x in chains])
                mdl = self._chain_mapping_mdl.Select(query)

                # chem mapping which is reduced to the respective chains
                chem_mapping = list()
                for m in self._chem_mapping:
                    chem_mapping.append([x for x in m if x in chains]) 

                self._get_repr_input[mdl_ligand.handle.hash_code] = \
                (mdl, chem_mapping)

            else:
                self._get_repr_input[mdl_ligand.handle.hash_code] = \
                (self._chain_mapping_mdl.CreateEmptyView(),
                 [list() for _ in self._chem_mapping])

        return (self._get_repr_input[mdl_ligand.hash_code][1],
                self._chem_group_alns,
                self._mdl_chains_without_chem_mapping,
                self._get_repr_input[mdl_ligand.hash_code][0])


def SCRMSD(model_ligand, target_ligand, transformation=geom.Mat4(),
           substructure_match=False, max_symmetries=1e6):
    """Calculate symmetry-corrected RMSD.

    Binding site superposition must be computed separately and passed as
    `transformation`.

    :param model_ligand: The model ligand
    :type model_ligand: :class:`ost.mol.ResidueHandle` or
                        :class:`ost.mol.ResidueView`
    :param target_ligand: The target ligand
    :type target_ligand: :class:`ost.mol.ResidueHandle` or
                         :class:`ost.mol.ResidueView`
    :param transformation: Optional transformation to apply on each atom
                           position of model_ligand.
    :type transformation: :class:`ost.geom.Mat4`
    :param substructure_match: Set this to True to allow partial target
                               ligand.
    :type substructure_match: :class:`bool`
    :param max_symmetries: If more than that many isomorphisms exist, raise
      a :class:`TooManySymmetriesError`. This can only be assessed by
      generating at least that many isomorphisms and can take some time.
    :type max_symmetries: :class:`int`
    :rtype: :class:`float`
    :raises: :class:`ost.mol.alg.ligand_scoring_base.NoSymmetryError` when no
             symmetry can be found,
             :class:`ost.mol.alg.ligand_scoring_base.DisconnectedGraphError`
             when ligand graph is disconnected,
             :class:`ost.mol.alg.ligand_scoring_base.TooManySymmetriesError`
             when more than *max_symmetries* isomorphisms are found.
    """

    symmetries = ligand_scoring_base.ComputeSymmetries(model_ligand,
                                        target_ligand,
                                        substructure_match=substructure_match,
                                        by_atom_index=True,
                                        max_symmetries=max_symmetries)
    return _SCRMSD_symmetries(symmetries, model_ligand, target_ligand,
                              transformation)


def _SCRMSD_symmetries(symmetries, model_ligand, target_ligand, 
                       transformation):
    """Compute SCRMSD with pre-computed symmetries. Internal. """

    # setup numpy positions for model ligand and apply transformation
    mdl_ligand_pos = np.ones((model_ligand.GetAtomCount(), 4))
    for a_idx, a in enumerate(model_ligand.atoms):
        p = a.GetPos()
        mdl_ligand_pos[a_idx, 0] = p[0]
        mdl_ligand_pos[a_idx, 1] = p[1]
        mdl_ligand_pos[a_idx, 2] = p[2]
    np_transformation = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            np_transformation[i,j] = transformation[i,j]
    mdl_ligand_pos = mdl_ligand_pos.dot(np_transformation.T)[:,:3]

    # setup numpy positions for target ligand
    trg_ligand_pos = np.zeros((target_ligand.GetAtomCount(), 3))
    for a_idx, a in enumerate(target_ligand.atoms):
        p = a.GetPos()
        trg_ligand_pos[a_idx, 0] = p[0]
        trg_ligand_pos[a_idx, 1] = p[1]
        trg_ligand_pos[a_idx, 2] = p[2]

    # position matrices to iterate symmetries
    # there is a guarantee that
    # target_ligand.GetAtomCount() <= model_ligand.GetAtomCount()
    # and that each target ligand atom is part of every symmetry
    # => target_ligand.GetAtomCount() is size of both position matrices
    rmsd_mdl_pos = np.zeros((target_ligand.GetAtomCount(), 3))
    rmsd_trg_pos = np.zeros((target_ligand.GetAtomCount(), 3))

    # iterate symmetries and find the one with lowest RMSD
    best_rmsd = np.inf
    for i, (trg_sym, mdl_sym) in enumerate(symmetries):
        for idx, (mdl_anum, trg_anum) in enumerate(zip(mdl_sym, trg_sym)):
            rmsd_mdl_pos[idx,:] = mdl_ligand_pos[mdl_anum, :]
            rmsd_trg_pos[idx,:] = trg_ligand_pos[trg_anum, :]
        rmsd = np.sqrt(((rmsd_mdl_pos - rmsd_trg_pos)**2).sum(-1).mean())
        if rmsd < best_rmsd:
            best_rmsd = rmsd

    return best_rmsd

# specify public interface
__all__ = ('SCRMSDScorer', 'SCRMSD')
