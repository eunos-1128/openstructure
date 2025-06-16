import numpy as np

from ost import LogWarning, LogInfo
from ost import geom
from ost import mol
from ost import seq

from ost.mol.alg import lddt
from ost.mol.alg import chain_mapping
from ost.mol.alg import ligand_scoring_base

class LDDTPLIScorer(ligand_scoring_base.LigandScorer):
    """ :class:`LigandScorer` implementing LDDT-PLI.

    LDDT-PLI is an LDDT score considering contacts between ligand and
    receptor. Where receptor consists of protein and nucleic acid chains that
    pass the criteria for :class:`chain mapping <ost.mol.alg.chain_mapping>`.
    This means ignoring other ligands, waters, short polymers as well as any
    incorrectly connected chains that may be in proximity.

    :class:`LDDTPLIScorer` computes a score for a specific pair of target/model
    ligands. Given a target/model ligand pair, all possible mappings of
    model chains onto their chemically equivalent target chains are enumerated.
    For each of these enumerations, all possible symmetries, i.e. atom-atom
    assignments of the ligand as given by :class:`LigandScorer`, are evaluated
    and an LDDT-PLI score is computed. The best possible LDDT-PLI score is
    returned.

    The LDDT-PLI score is a variant of LDDT with a custom inclusion radius
    (`lddt_pli_radius`), no stereochemistry checks, and which penalizes
    contacts added in the model within `lddt_pli_radius` by default
    (can be changed with the `add_mdl_contacts` flag) but only if the involved
    atoms can be mapped to the target. This is a requirement to
    1) extract the respective reference distance from the target
    2) avoid usage of contacts for which we have no experimental evidence.
    One special case are contacts from chains that are not mapped to the target
    binding site. It is very well possible that we have experimental evidence
    for this chain though its just too far away from the target binding site.
    We therefore try to map these contacts to the chain in the target with
    equivalent sequence that is closest to the target binding site. If the
    respective atoms can be mapped there, the contact is considered not
    fulfilled and added as penalty.

    Populates :attr:`LigandScorer.aux_data` with following :class:`dict` keys:

    * lddt_pli: The LDDT-PLI score
    * lddt_pli_n_contacts: Number of contacts considered in LDDT computation
    * target_ligand: The actual target ligand for which the score was computed
    * model_ligand: The actual model ligand for which the score was computed
    * chain_mapping: :class:`dict` with a chain mapping of chains involved in
      binding site - key: trg chain name, value: mdl chain name

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
    :param lddt_pli_radius: LDDT inclusion radius for LDDT-PLI.
    :type lddt_pli_radius: :class:`float`
    :param add_mdl_contacts: Whether to penalize added model contacts.
    :type add_mdl_contacts: :class:`bool`
    :param lddt_pli_thresholds: Distance difference thresholds for LDDT.
    :type lddt_pli_thresholds: :class:`list` of :class:`float`
    :param lddt_pli_binding_site_radius: Pro param - dont use. Providing a value
                                         Restores behaviour from previous
                                         implementation that first extracted a
                                         binding site with strict distance
                                         threshold and computed LDDT-PLI only on
                                         those target residues whereas the
                                         current implementation includes every
                                         atom within *lddt_pli_radius*.
    :type lddt_pli_binding_site_radius: :class:`float`
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
    """

    def __init__(self, model, target, model_ligands, target_ligands,
                 resnum_alignments=False, rename_ligand_chain=False,
                 substructure_match=False, coverage_delta=0.2,
                 max_symmetries=1e4, lddt_pli_radius=6.0,
                 add_mdl_contacts=True,
                 lddt_pli_thresholds = [0.5, 1.0, 2.0, 4.0],
                 lddt_pli_binding_site_radius=None,
                 min_pep_length = 6,
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

        self.lddt_pli_radius = lddt_pli_radius
        self.add_mdl_contacts = add_mdl_contacts
        self.lddt_pli_thresholds = lddt_pli_thresholds
        self.lddt_pli_binding_site_radius = lddt_pli_binding_site_radius

        # lazily precomputed variables to speedup lddt-pli computation
        self._lddt_pli_target_data = dict()
        self._lddt_pli_model_data = dict()
        self.__mappable_atoms = None

        # update state decoding from parent with subclass specific stuff
        self.state_decoding[10] = ("no_contact",
                                   "There were no LDDT contacts between the "
                                   "binding site and the ligand, and LDDT-PLI "
                                   "is undefined.")
        self.state_decoding[20] = ("unknown",
                                   "Unknown error occured in LDDTPLIScorer")

    def _compute(self, symmetries, target_ligand, model_ligand):
        """ Implements interface from parent
        """
        if self.add_mdl_contacts:
            LogInfo("Computing LDDT-PLI with added model contacts")
            result = self._compute_lddt_pli_add_mdl_contacts(symmetries,
                                                             target_ligand,
                                                             model_ligand)
        else:
            LogInfo("Computing LDDT-PLI without added model contacts")
            result = self._compute_lddt_pli_classic(symmetries,
                                                    target_ligand,
                                                    model_ligand)

        pair_state = 0
        score = result["lddt_pli"]

        if score is None or np.isnan(score):
            if result["lddt_pli_n_contacts"] == 0:
                # it's a space ship!
                pair_state = 10
            else:
                # unknwon error state
                pair_state = 20

        # the ligands get a zero-state...
        target_ligand_state = 0
        model_ligand_state = 0

        return (score, pair_state, target_ligand_state, model_ligand_state,
                result)

    def _score_dir(self):
        """ Implements interface from parent
        """
        return '+'

    def _compute_lddt_pli_add_mdl_contacts(self, symmetries, target_ligand,
                                           model_ligand):

        ###############################
        # Get stuff from model/target #
        ###############################

        trg_residues, trg_bs, trg_chains, trg_ligand_chain, \
        trg_ligand_res, scorer, chem_groups = \
        self._lddt_pli_get_trg_data(target_ligand)

        trg_bs_center = trg_bs.geometric_center

        # Copy to make sure that we don't change anything on underlying
        # references
        # This is not strictly necessary in the current implementation but
        # hey, maybe it avoids hard to debug errors when someone changes things
        ref_indices = [a.copy() for a in scorer.ref_indices_ic]
        ref_distances = [a.copy() for a in scorer.ref_distances_ic]

        # distance hacking... remove any interchain distance except the ones
        # with the ligand
        ligand_start_idx = scorer.chain_start_indices[-1]
        for at_idx in range(ligand_start_idx):
            mask = ref_indices[at_idx] >= ligand_start_idx
            ref_indices[at_idx] = ref_indices[at_idx][mask]
            ref_distances[at_idx] = ref_distances[at_idx][mask]

        mdl_residues, mdl_bs, mdl_chains, mdl_ligand_chain, mdl_ligand_res, \
        chem_mapping = self._lddt_pli_get_mdl_data(model_ligand)

        ####################
        # Setup alignments #
        ####################

        # ref_mdl_alns refers to full chain mapper trg and mdl structures
        # => need to adapt mdl sequence that only contain residues in contact
        #    with ligand
        cut_ref_mdl_alns = self._lddt_pli_cut_ref_mdl_alns(chem_groups,
                                                           chem_mapping,
                                                           mdl_bs, trg_bs)

        ########################################
        # Setup cache for added model contacts #
        ########################################

        # get each chain mapping that we ever observe in scoring
        chain_mappings = list(chain_mapping._ChainMappings(chem_groups,
                                                           chem_mapping))

        # for each mdl ligand atom, we collect all trg ligand atoms that are
        # ever mapped onto it given *symmetries*
        ligand_atom_mappings = [set() for a in mdl_ligand_res.atoms]
        for (trg_sym, mdl_sym) in symmetries:
            for trg_i, mdl_i in zip(trg_sym, mdl_sym):
                ligand_atom_mappings[mdl_i].add(trg_i)

        mdl_ligand_pos = np.zeros((mdl_ligand_res.GetAtomCount(), 3))
        for a_idx, a in enumerate(mdl_ligand_res.atoms):
            p = a.GetPos()
            mdl_ligand_pos[a_idx, 0] = p[0]
            mdl_ligand_pos[a_idx, 1] = p[1]
            mdl_ligand_pos[a_idx, 2] = p[2]

        trg_ligand_pos = np.zeros((trg_ligand_res.GetAtomCount(), 3))
        for a_idx, a in enumerate(trg_ligand_res.atoms):
            p = a.GetPos()
            trg_ligand_pos[a_idx, 0] = p[0]
            trg_ligand_pos[a_idx, 1] = p[1]
            trg_ligand_pos[a_idx, 2] = p[2]

        mdl_lig_hashes = [a.hash_code for a in mdl_ligand_res.atoms]

        symmetric_atoms = np.asarray(sorted(list(scorer.symmetric_atoms)),
                                     dtype=np.int64)

        # two caches to cache things for each chain mapping => lists
        # of len(chain_mappings)
        # 
        # In principle we're caching for each trg/mdl ligand atom pair all
        # information to update ref_indices/ref_distances and resolving the
        # symmetries of the binding site.
        # in detail: each list entry in *scoring_cache* is a dict with
        # key: (mdl_lig_at_idx, trg_lig_at_idx)
        # value: tuple with 4 elements - 1: indices of atoms representing added
        # contacts relative to overall inexing scheme in scorer 2: the
        # respective distances 3: the same but only containing indices towards
        # atoms of the binding site that are considered symmetric 4: the
        # respective indices.
        # each list entry in *penalty_cache* is a list of len N mdl lig atoms.
        # For each mdl lig at it contains a penalty for this mdl lig at. That
        # means the number of contacts in the mdl binding site that can
        # directly be mapped to the target given the local chain mapping but
        # are not present in the target binding site, i.e. interacting atoms are
        # too far away.
        scoring_cache = list()
        penalty_cache = list()

        for mapping in chain_mappings:

            # flat mapping with mdl chain names as key
            flat_mapping = dict()
            for trg_chem_group, mdl_chem_group in zip(chem_groups, mapping):
                for a,b in zip(trg_chem_group, mdl_chem_group):
                    if a is not None and b is not None:
                        flat_mapping[b] = a

            # for each mdl bs atom (as atom hash), the trg bs atoms (as index in
            # scorer)
            bs_atom_mapping = dict()
            for mdl_cname, ref_cname in flat_mapping.items():
                aln = cut_ref_mdl_alns[(ref_cname, mdl_cname)]
                ref_ch = trg_bs.Select(f"cname={mol.QueryQuoteName(ref_cname)}")
                mdl_ch = mdl_bs.Select(f"cname={mol.QueryQuoteName(mdl_cname)}")
                aln.AttachView(0, ref_ch)
                aln.AttachView(1, mdl_ch)
                for col in aln:
                    ref_r = col.GetResidue(0)
                    mdl_r = col.GetResidue(1)
                    if ref_r.IsValid() and mdl_r.IsValid():
                        for mdl_a in mdl_r.atoms:
                            ref_a = ref_r.FindAtom(mdl_a.GetName())
                            if ref_a.IsValid():
                                ref_h = ref_a.handle.hash_code
                                if ref_h in scorer.atom_indices:
                                    mdl_h = mdl_a.handle.hash_code
                                    bs_atom_mapping[mdl_h] = \
                                    scorer.atom_indices[ref_h]

            cache = dict()
            n_penalties = list()

            for mdl_a_idx, mdl_a in enumerate(mdl_ligand_res.atoms):
                n_penalty = 0
                trg_bs_indices = list()
                close_a = mdl_bs.FindWithin(mdl_a.GetPos(),
                                            self.lddt_pli_radius)
                for a in close_a:
                    mdl_a_hash_code = a.hash_code
                    if mdl_a_hash_code in bs_atom_mapping:
                        trg_bs_indices.append(bs_atom_mapping[mdl_a_hash_code])
                    elif mdl_a_hash_code not in mdl_lig_hashes:
                        if a.GetChain().GetName() in flat_mapping:
                            # Its in a mapped chain
                            at_key = (a.GetResidue().GetNumber(), a.name)
                            cname = a.GetChain().name
                            cname_key = (flat_mapping[cname], cname)
                            if at_key in self._mappable_atoms[cname_key]:
                                # Its a contact in the model but not part of
                                # trg_bs. It can still be mapped using the
                                # global mdl_ch/ref_ch alignment
                                # d in ref > self.lddt_pli_radius + max_thresh
                                # => guaranteed to be non-fulfilled contact
                                n_penalty += 1

                n_penalties.append(n_penalty)

                trg_bs_indices = np.asarray(sorted(trg_bs_indices))

                for trg_a_idx in ligand_atom_mappings[mdl_a_idx]:
                    # mask selects entries in trg_bs_indices that are not yet
                    # part of classic LDDT ref_indices for atom at trg_a_idx
                    # => added mdl contacts
                    mask = np.isin(trg_bs_indices,
                                   ref_indices[ligand_start_idx + trg_a_idx],
                                   assume_unique=True, invert=True)
                    added_indices = np.asarray([], dtype=np.int64)
                    added_distances = np.asarray([], dtype=np.float64)
                    if np.sum(mask) > 0:
                        # compute ref distances on reference positions
                        added_indices = trg_bs_indices[mask]
                        tmp = scorer.positions.take(added_indices, axis=0)
                        np.subtract(tmp, trg_ligand_pos[trg_a_idx][None, :],
                                    out=tmp)
                        np.square(tmp, out=tmp)
                        tmp = tmp.sum(axis=1)
                        np.sqrt(tmp, out=tmp)
                        added_distances = tmp

                    # extract the distances towards bs atoms that are symmetric
                    sym_mask = np.isin(added_indices, symmetric_atoms,
                                       assume_unique=True)

                    cache[(mdl_a_idx, trg_a_idx)] = (added_indices,
                                                     added_distances,
                                                     added_indices[sym_mask],
                                                     added_distances[sym_mask])

            scoring_cache.append(cache)
            penalty_cache.append(n_penalties)

        # cache for model contacts towards non mapped trg chains - this is
        # relevant for self._lddt_pli_unmapped_chain_penalty
        # key: tuple in form (trg_ch, mdl_ch)
        # value: yet another dict with
        #   key: ligand_atom_hash
        #   value: n contacts towards respective trg chain that can be mapped
        non_mapped_cache = dict()

        ###############################################################
        # compute LDDT for all possible chain mappings and symmetries #
        ###############################################################

        best_score = -1.0
        best_result = {"lddt_pli": None,
                       "lddt_pli_n_contacts": 0,
                       "chain_mapping": None}

        # dummy alignment for ligand chains which is needed as input later on
        ligand_aln = seq.CreateAlignment()
        trg_s = seq.CreateSequence(trg_ligand_chain.name,
                                   trg_ligand_res.GetOneLetterCode())
        mdl_s = seq.CreateSequence(mdl_ligand_chain.name,
                                   mdl_ligand_res.GetOneLetterCode())
        ligand_aln.AddSequence(trg_s)
        ligand_aln.AddSequence(mdl_s)
        ligand_at_indices =  list(range(ligand_start_idx, scorer.n_atoms))

        sym_idx_collector = [None] * scorer.n_atoms
        sym_dist_collector = [None] * scorer.n_atoms

        for mapping, s_cache, p_cache in zip(chain_mappings, scoring_cache,
                                             penalty_cache):

            lddt_chain_mapping = dict()
            lddt_alns = dict()
            for ref_chem_group, mdl_chem_group in zip(chem_groups, mapping):
                for ref_ch, mdl_ch in zip(ref_chem_group, mdl_chem_group):
                    # some mdl chains can be None
                    if mdl_ch is not None:
                        lddt_chain_mapping[mdl_ch] = ref_ch
                        lddt_alns[mdl_ch] = cut_ref_mdl_alns[(ref_ch, mdl_ch)]

            # add ligand to lddt_chain_mapping/lddt_alns
            lddt_chain_mapping[mdl_ligand_chain.name] = trg_ligand_chain.name
            lddt_alns[mdl_ligand_chain.name] = ligand_aln

            # already process model, positions will be manually hacked for each
            # symmetry - small overhead for variables that are thrown away here
            pos, _, _, _, _, _, lddt_symmetries = \
            scorer._ProcessModel(mdl_bs, lddt_chain_mapping,
                                 residue_mapping = lddt_alns,
                                 nirvana_dist = self.lddt_pli_radius + max(self.lddt_pli_thresholds),
                                 check_resnames = False)

            # estimate a penalty for unsatisfied model contacts from chains
            # that are not in the local trg binding site, but can be mapped in
            # the target.
            # We're using the trg chain with the closest geometric center to
            # the trg binding site that can be mapped to the mdl chain
            # according the chem mapping. An alternative would be to search for
            # the target chain with the minimal number of additional contacts.
            # There is not good solution for this problem...
            unmapped_chains = list()
            already_mapped = set()
            for mdl_ch in mdl_chains:
                if mdl_ch not in lddt_chain_mapping:

                    if mdl_ch in self._mdl_chains_without_chem_mapping:
                        # this mdl chain does not map to any trg chain
                        continue

                    # check which chain in trg is closest
                    chem_grp_idx = None
                    for i, m in enumerate(self._chem_mapping):
                        if mdl_ch in m:
                            chem_grp_idx = i
                            break
                    if chem_grp_idx is None:
                        raise RuntimeError("This should never happen... "
                                           "ask Gabriel...")
                    closest_ch = None
                    closest_dist = None
                    for trg_ch in self._chain_mapper.chem_groups[chem_grp_idx]:
                        if trg_ch not in lddt_chain_mapping.values():
                            if trg_ch not in already_mapped:
                                ch = self._chain_mapper.target.FindChain(trg_ch) 
                                c = ch.geometric_center
                                d = geom.Distance(trg_bs_center, c)
                                if closest_dist is None or d < closest_dist:
                                    closest_dist = d
                                    closest_ch = trg_ch
                    if closest_ch is not None:
                        unmapped_chains.append((closest_ch, mdl_ch))
                        already_mapped.add(closest_ch)

            for (trg_sym, mdl_sym) in symmetries:

                # update positions
                for mdl_i, trg_i in zip(mdl_sym, trg_sym):
                    pos[ligand_start_idx + trg_i, :] = mdl_ligand_pos[mdl_i, :]

                # start new ref_indices/ref_distances from original values
                funky_ref_indices = [np.copy(a) for a in ref_indices]
                funky_ref_distances = [np.copy(a) for a in ref_distances]

                # The only distances from the binding site towards the ligand
                # we care about are the ones from the symmetric atoms to
                # correctly compute scorer._ResolveSymmetries.
                # We collect them while updating distances from added mdl
                # contacts
                for idx in symmetric_atoms:
                    sym_idx_collector[idx] = list()
                    sym_dist_collector[idx] = list()

                # add data from added mdl contacts cache
                added_penalty = 0
                for mdl_i, trg_i in zip(mdl_sym, trg_sym):
                    added_penalty += p_cache[mdl_i]
                    cache = s_cache[mdl_i, trg_i]
                    full_trg_i = ligand_start_idx + trg_i
                    funky_ref_indices[full_trg_i] = \
                    np.append(funky_ref_indices[full_trg_i], cache[0])
                    funky_ref_distances[full_trg_i] = \
                    np.append(funky_ref_distances[full_trg_i], cache[1])
                    for idx, d in zip(cache[2], cache[3]):
                        sym_idx_collector[idx].append(full_trg_i)
                        sym_dist_collector[idx].append(d)

                for idx in symmetric_atoms:
                    funky_ref_indices[idx] = \
                    np.append(funky_ref_indices[idx],
                              np.asarray(sym_idx_collector[idx],
                                         dtype=np.int64))
                    funky_ref_distances[idx] = \
                    np.append(funky_ref_distances[idx],
                              np.asarray(sym_dist_collector[idx],
                                         dtype=np.float64))

                # we can pass funky_ref_indices/funky_ref_distances as
                # sym_ref_indices/sym_ref_distances in
                # scorer._ResolveSymmetries as we only have distances of the bs
                # to the ligand and ligand atoms are "non-symmetric"
                scorer._ResolveSymmetries(pos, self.lddt_pli_thresholds,
                                          lddt_symmetries,
                                          funky_ref_indices,
                                          funky_ref_distances)

                N = sum([len(funky_ref_indices[i]) for i in ligand_at_indices])
                N += added_penalty

                # collect number of expected contacts which can be mapped
                if len(unmapped_chains) > 0:
                    N += self._lddt_pli_unmapped_chain_penalty(unmapped_chains,
                                                               non_mapped_cache,
                                                               mdl_bs,
                                                               mdl_ligand_res,
                                                               mdl_sym)

                conserved = np.sum(scorer._EvalAtoms(pos, ligand_at_indices,
                                                     self.lddt_pli_thresholds,
                                                     funky_ref_indices,
                                                     funky_ref_distances),
                                   axis=0)
                score = None
                if N > 0:
                    score = np.mean(conserved/N)

                if score is not None and score > best_score:
                    best_score = score
                    save_chain_mapping = dict(lddt_chain_mapping)
                    del save_chain_mapping[mdl_ligand_chain.name]
                    best_result = {"lddt_pli": score,
                                   "lddt_pli_n_contacts": N,
                                   "chain_mapping": save_chain_mapping}

        # fill misc info to result object
        best_result["target_ligand"] = target_ligand
        best_result["model_ligand"] = model_ligand

        return best_result


    def _compute_lddt_pli_classic(self, symmetries, target_ligand,
                                  model_ligand):

        ###############################
        # Get stuff from model/target #
        ###############################

        max_r = None
        if self.lddt_pli_binding_site_radius:
            max_r = self.lddt_pli_binding_site_radius

        trg_residues, trg_bs, trg_chains, trg_ligand_chain, \
        trg_ligand_res, scorer, chem_groups = \
        self._lddt_pli_get_trg_data(target_ligand, max_r = max_r)

        # Copy to make sure that we don't change anything on underlying
        # references
        # This is not strictly necessary in the current implementation but
        # hey, maybe it avoids hard to debug errors when someone changes things
        ref_indices = [a.copy() for a in scorer.ref_indices_ic]
        ref_distances = [a.copy() for a in scorer.ref_distances_ic]

        # no matter what mapping/symmetries, the number of expected
        # contacts stays the same
        ligand_start_idx = scorer.chain_start_indices[-1]
        ligand_at_indices = list(range(ligand_start_idx, scorer.n_atoms))
        n_exp = sum([len(ref_indices[i]) for i in ligand_at_indices])

        mdl_residues, mdl_bs, mdl_chains, mdl_ligand_chain, mdl_ligand_res, \
        chem_mapping = self._lddt_pli_get_mdl_data(model_ligand)

        if n_exp == 0:
            # no contacts... nothing to compute...
            return {"lddt_pli": None,
                    "lddt_pli_n_contacts": 0,
                    "chain_mapping": None,
                    "target_ligand": target_ligand,
                    "model_ligand": model_ligand}

        # Distance hacking... remove any interchain distance except the ones
        # with the ligand
        for at_idx in range(ligand_start_idx):
            mask = ref_indices[at_idx] >= ligand_start_idx
            ref_indices[at_idx] = ref_indices[at_idx][mask]
            ref_distances[at_idx] = ref_distances[at_idx][mask]

        ####################
        # Setup alignments #
        ####################

        # ref_mdl_alns refers to full chain mapper trg and mdl structures
        # => need to adapt mdl sequence that only contain residues in contact
        #    with ligand
        cut_ref_mdl_alns = self._lddt_pli_cut_ref_mdl_alns(chem_groups,
                                                           chem_mapping,
                                                           mdl_bs, trg_bs)

        ###############################################################
        # compute LDDT for all possible chain mappings and symmetries #
        ###############################################################

        best_score = -1.0

        # dummy alignment for ligand chains which is needed as input later on
        l_aln = seq.CreateAlignment()
        l_aln.AddSequence(seq.CreateSequence(trg_ligand_chain.name,
                                             trg_ligand_res.GetOneLetterCode()))
        l_aln.AddSequence(seq.CreateSequence(mdl_ligand_chain.name,
                                             mdl_ligand_res.GetOneLetterCode()))

        mdl_ligand_pos = np.zeros((model_ligand.GetAtomCount(), 3))
        for a_idx, a in enumerate(model_ligand.atoms):
            p = a.GetPos()
            mdl_ligand_pos[a_idx, 0] = p[0]
            mdl_ligand_pos[a_idx, 1] = p[1]
            mdl_ligand_pos[a_idx, 2] = p[2]

        for mapping in chain_mapping._ChainMappings(chem_groups, chem_mapping):

            lddt_chain_mapping = dict()
            lddt_alns = dict()
            for ref_chem_group, mdl_chem_group in zip(chem_groups, mapping):
                for ref_ch, mdl_ch in zip(ref_chem_group, mdl_chem_group):
                    # some mdl chains can be None
                    if mdl_ch is not None:
                        lddt_chain_mapping[mdl_ch] = ref_ch
                        lddt_alns[mdl_ch] = cut_ref_mdl_alns[(ref_ch, mdl_ch)]

            # add ligand to lddt_chain_mapping/lddt_alns
            lddt_chain_mapping[mdl_ligand_chain.name] = trg_ligand_chain.name
            lddt_alns[mdl_ligand_chain.name] = l_aln

            # already process model, positions will be manually hacked for each
            # symmetry - small overhead for variables that are thrown away here
            pos, _, _, _, _, _, lddt_symmetries = \
            scorer._ProcessModel(mdl_bs, lddt_chain_mapping,
                                 residue_mapping = lddt_alns,
                                 nirvana_dist = self.lddt_pli_radius + max(self.lddt_pli_thresholds),
                                 check_resnames = False)

            for (trg_sym, mdl_sym) in symmetries:
                for mdl_i, trg_i in zip(mdl_sym, trg_sym):
                    pos[ligand_start_idx + trg_i, :] = mdl_ligand_pos[mdl_i, :]
                # we can pass ref_indices/ref_distances as
                # sym_ref_indices/sym_ref_distances in
                # scorer._ResolveSymmetries as we only have distances of the bs
                # to the ligand and ligand atoms are "non-symmetric"
                scorer._ResolveSymmetries(pos, self.lddt_pli_thresholds,
                                          lddt_symmetries,
                                          ref_indices,
                                          ref_distances)
                # compute number of conserved distances for ligand atoms
                conserved = np.sum(scorer._EvalAtoms(pos, ligand_at_indices,
                                                     self.lddt_pli_thresholds,
                                                     ref_indices,
                                                     ref_distances), axis=0)
                score = np.mean(conserved/n_exp)

                if score > best_score:
                    best_score = score
                    save_chain_mapping = dict(lddt_chain_mapping)
                    del save_chain_mapping[mdl_ligand_chain.name]
                    best_result = {"lddt_pli": score,
                                   "chain_mapping": save_chain_mapping}

        # fill misc info to result object
        best_result["lddt_pli_n_contacts"] = n_exp
        best_result["target_ligand"] = target_ligand
        best_result["model_ligand"] = model_ligand

        return best_result

    def _lddt_pli_unmapped_chain_penalty(self, unmapped_chains,
                                         non_mapped_cache,
                                         mdl_bs,
                                         mdl_ligand_res,
                                         mdl_sym):

        n_exp = 0
        for ch_tuple in unmapped_chains:
            if ch_tuple not in non_mapped_cache:
                # for each ligand atom, we count the number of mappable atoms
                # within lddt_pli_radius
                counts = dict()
                # the select statement also excludes the ligand in mdl_bs
                # as it resides in a separate chain
                mdl_cname = ch_tuple[1]
                query = "cname=" + mol.QueryQuoteName(mdl_cname)
                mdl_bs_ch = mdl_bs.Select(query)
                for a in mdl_ligand_res.atoms:
                    close_atoms = \
                    mdl_bs_ch.FindWithin(a.GetPos(), self.lddt_pli_radius)
                    N = 0
                    for close_a in close_atoms:
                        at_key = (close_a.GetResidue().GetNumber(),
                                  close_a.GetName())
                        if at_key in self._mappable_atoms[ch_tuple]:
                            N += 1
                    counts[a.hash_code] = N

                # fill cache
                non_mapped_cache[ch_tuple] = counts

            # add number of mdl contacts which can be mapped to target
            # as non-fulfilled contacts
            counts = non_mapped_cache[ch_tuple]
            lig_hash_codes = [a.hash_code for a in mdl_ligand_res.atoms]
            for i in mdl_sym:
                n_exp += counts[lig_hash_codes[i]]

        return n_exp


    def _lddt_pli_get_mdl_data(self, model_ligand):
        if model_ligand not in self._lddt_pli_model_data:

            mdl = self._chain_mapping_mdl

            mdl_residues = set()
            for at in model_ligand.atoms:
                close_atoms = mdl.FindWithin(at.GetPos(), self.lddt_pli_radius)
                for close_at in close_atoms:
                    mdl_residues.add(close_at.GetResidue())

            max_r = self.lddt_pli_radius + max(self.lddt_pli_thresholds)
            for r in mdl.residues:
                r.SetIntProp("bs", 0)
            for at in model_ligand.atoms:
                close_atoms = mdl.FindWithin(at.GetPos(), max_r)
                for close_at in close_atoms:
                    close_at.GetResidue().SetIntProp("bs", 1)

            mdl_bs = mol.CreateEntityFromView(mdl.Select("grbs:0=1"), True)
            mdl_chains = set([ch.name for ch in mdl_bs.chains])

            mdl_editor = mdl_bs.EditXCS(mol.BUFFERED_EDIT)
            mdl_ligand_chain = None
            for cname in ["hugo_the_cat_terminator", "ida_the_cheese_monster"]:
                try:
                    # I'm pretty sure, one of these chain names is not there...
                    mdl_ligand_chain = mdl_editor.InsertChain(cname)
                    break
                except:
                    pass
            if mdl_ligand_chain is None:
                raise RuntimeError("Fuck this, I'm out...")
            mdl_ligand_res = mdl_editor.AppendResidue(mdl_ligand_chain,
                                                      model_ligand,
                                                      deep=True)
            mdl_editor.RenameResidue(mdl_ligand_res, "LIG")
            mdl_editor.SetResidueNumber(mdl_ligand_res, mol.ResNum(1))

            chem_mapping = list()
            for m in self._chem_mapping:
                chem_mapping.append([x for x in m if x in mdl_chains])

            self._lddt_pli_model_data[model_ligand] = (mdl_residues,
                                                       mdl_bs,
                                                       mdl_chains,
                                                       mdl_ligand_chain,
                                                       mdl_ligand_res,
                                                       chem_mapping)

        return self._lddt_pli_model_data[model_ligand]


    def _lddt_pli_get_trg_data(self, target_ligand, max_r = None):
        if target_ligand not in self._lddt_pli_target_data:

            trg = self._chain_mapper.target

            if max_r is None:
                max_r = self.lddt_pli_radius + max(self.lddt_pli_thresholds)

            trg_residues = set()
            for at in target_ligand.atoms:
                close_atoms = trg.FindWithin(at.GetPos(), max_r)
                for close_at in close_atoms:
                    trg_residues.add(close_at.GetResidue())

            for r in trg.residues:
                r.SetIntProp("bs", 0)

            for r in trg_residues:
                r.SetIntProp("bs", 1)

            trg_bs = mol.CreateEntityFromView(trg.Select("grbs:0=1"), True)
            trg_chains = set([ch.name for ch in trg_bs.chains])

            trg_editor = trg_bs.EditXCS(mol.BUFFERED_EDIT)
            trg_ligand_chain = None
            for cname in ["hugo_the_cat_terminator", "ida_the_cheese_monster"]:
                try:
                    # I'm pretty sure, one of these chain names is not there yet
                    trg_ligand_chain = trg_editor.InsertChain(cname)
                    break
                except:
                    pass
            if trg_ligand_chain is None:
                raise RuntimeError("Fuck this, I'm out...")

            trg_ligand_res = trg_editor.AppendResidue(trg_ligand_chain,
                                                      target_ligand,
                                                      deep=True)
            trg_editor.RenameResidue(trg_ligand_res, "LIG")
            trg_editor.SetResidueNumber(trg_ligand_res, mol.ResNum(1))

            compound_name = trg_ligand_res.name
            compound = lddt.CustomCompound.FromResidue(trg_ligand_res)
            custom_compounds = {compound_name: compound}

            scorer = lddt.lDDTScorer(trg_bs,
                                     custom_compounds = custom_compounds,
                                     inclusion_radius = self.lddt_pli_radius)

            chem_groups = list()
            for g in self._chain_mapper.chem_groups:
                chem_groups.append([x for x in g if x in trg_chains])

            self._lddt_pli_target_data[target_ligand] = (trg_residues,
                                                         trg_bs,
                                                         trg_chains,
                                                         trg_ligand_chain,
                                                         trg_ligand_res,
                                                         scorer,
                                                         chem_groups)

        return self._lddt_pli_target_data[target_ligand]


    def _lddt_pli_cut_ref_mdl_alns(self, chem_groups, chem_mapping, mdl_bs,
                                   ref_bs):
        cut_ref_mdl_alns = dict()
        for ref_chem_group, mdl_chem_group in zip(chem_groups, chem_mapping):
            for ref_ch in ref_chem_group:

                ref_bs_chain = ref_bs.FindChain(ref_ch)
                query = "cname=" + mol.QueryQuoteName(ref_ch)
                ref_view = self._chain_mapper.target.Select(query)

                for mdl_ch in mdl_chem_group:
                    aln = self._ref_mdl_alns[(ref_ch, mdl_ch)]

                    aln.AttachView(0, ref_view)

                    mdl_bs_chain = mdl_bs.FindChain(mdl_ch)
                    query = "cname=" + mol.QueryQuoteName(mdl_ch)
                    aln.AttachView(1, self._chain_mapping_mdl.Select(query))

                    cut_mdl_seq = ['-'] * aln.GetLength()
                    cut_ref_seq = ['-'] * aln.GetLength()
                    for i, col in enumerate(aln):

                       # check ref residue
                        r = col.GetResidue(0)
                        if r.IsValid():
                            bs_r = ref_bs_chain.FindResidue(r.GetNumber())
                            if bs_r.IsValid():
                                cut_ref_seq[i] = col[0]

                        # check mdl residue
                        r = col.GetResidue(1)
                        if r.IsValid():
                            bs_r = mdl_bs_chain.FindResidue(r.GetNumber())
                            if bs_r.IsValid():
                                cut_mdl_seq[i] = col[1]

                    cut_ref_seq = ''.join(cut_ref_seq)         
                    cut_mdl_seq = ''.join(cut_mdl_seq)         
                    cut_aln = seq.CreateAlignment()
                    cut_aln.AddSequence(seq.CreateSequence(ref_ch, cut_ref_seq))
                    cut_aln.AddSequence(seq.CreateSequence(mdl_ch, cut_mdl_seq))
                    cut_ref_mdl_alns[(ref_ch, mdl_ch)] = cut_aln
        return cut_ref_mdl_alns

    @property
    def _mappable_atoms(self):
        """ Stores mappable atoms given a chain mapping

        Store for each ref_ch,mdl_ch pair all mdl atoms that can be
        mapped. Don't store mappable atoms as hashes but rather as tuple
        (mdl_r.GetNumber(), mdl_a.GetName()). Reason for that is that one might
        operate on Copied EntityHandle objects without corresponding hashes.
        Given a tuple defining c_pair: (ref_cname, mdl_cname), one
        can check if a certain atom is mappable by evaluating:
        if (mdl_r.GetNumber(), mdl_a.GetName()) in self._mappable_atoms(c_pair)
        """
        if self.__mappable_atoms is None:
            self.__mappable_atoms = dict()
            for (ref_cname, mdl_cname), aln in self._ref_mdl_alns.items():
                self._mappable_atoms[(ref_cname, mdl_cname)] = set()
                ref_query = f"cname={mol.QueryQuoteName(ref_cname)}"
                mdl_query = f"cname={mol.QueryQuoteName(mdl_cname)}"
                ref_ch = self._chain_mapper.target.Select(ref_query)
                mdl_ch = self._chain_mapping_mdl.Select(mdl_query)
                aln.AttachView(0, ref_ch)
                aln.AttachView(1, mdl_ch)
                for col in aln:
                    ref_r = col.GetResidue(0)
                    mdl_r = col.GetResidue(1)
                    if ref_r.IsValid() and mdl_r.IsValid():
                        for mdl_a in mdl_r.atoms:
                            if ref_r.FindAtom(mdl_a.name).IsValid():
                                c_key = (ref_cname, mdl_cname)
                                at_key = (mdl_r.GetNumber(), mdl_a.name)
                                self.__mappable_atoms[c_key].add(at_key)

        return self.__mappable_atoms

# specify public interface
__all__ = ('LDDTPLIScorer',)
