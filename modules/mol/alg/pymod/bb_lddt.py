import itertools
import numpy as np
from scipy.spatial import distance

import time
from ost import mol

class BBlDDTEntity:
    """ Helper object for BBlDDT computation

    Holds structural information and getters for interacting chains, i.e.
    interfaces. Peptide residues are represented by their CA position
    and nucleotides by C3'.

    :param ent: Structure for BBlDDT score computation
    :type ent: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param contact_d: Pairwise distance of residues to be considered as contacts
    :type contact_d: :class:`float`
    """
    def __init__(self, ent, dist_thresh = 15.0,
                 dist_diff_thresholds = [0.5, 1.0, 2.0, 4.0]):
        pep_query = "(peptide=true and aname=\"CA\")"
        nuc_query = "(nucleotide=True and aname=\"C3'\")"
        self._view = ent.Select(" or ".join([pep_query, nuc_query]))
        self._dist_thresh = dist_thresh
        self._dist_diff_thresholds = dist_diff_thresholds

        # the following attributes will be lazily evaluated
        self._chain_names = None
        self._interacting_chains = None
        self._potentially_contributing_chains = None
        self._sequence = dict()
        self._pos = dict()
        self._pair_dist = dict()
        self._sc_dist = dict()
        self._n_pair_contacts = None
        self._n_sc_contacts = None
        self._n_contacts = None
        # min and max xyz for elements in pos used for fast collision
        # detection
        self._min_pos = dict()
        self._max_pos = dict()

    @property
    def view(self):
        """ Processed structure

        View that only contains representative atoms. That's CA for peptide
        residues and C3' for nucleotides.

        :type: :class:`ost.mol.EntityView`
        """
        return self._view

    @property
    def dist_thresh(self):
        """ Pairwise distance of residues to be considered as contacts

        Given at :class:`BBlDDTEntity` construction

        :type: :class:`float`
        """
        return self._dist_thresh

    @property
    def dist_diff_thresholds(self):
        """ Distance difference thresholds for lDDT computation

        Given at :class:`BBlDDTEntity` construction

        :type: :class:`list` of :class:`float`        
        """
        return self._dist_diff_thresholds

    @property
    def chain_names(self):
        """ Chain names in :attr:`~view`
 
        Names are sorted

        :type: :class:`list` of :class:`str`
        """
        if self._chain_names is None:
            self._chain_names = sorted([ch.name for ch in self.view.chains])
        return self._chain_names

    @property
    def interacting_chains(self):
        """ Pairs of chains in :attr:`~view` with at least one contact

        :type: :class:`list` of :class:`tuples`
        """
        if self._interacting_chains is None:
            # ugly hack: also computes self._n_pair_contacts
            # this assumption is made when computing the n_pair_contacts
            # attribute
            self._interacting_chains = list()
            self._n_pair_contacts = list()
            for x in itertools.combinations(self.chain_names, 2):
                if self.PotentialInteraction(x[0], x[1]):
                    n = np.count_nonzero(self.PairDist(x[0], x[1]) < self.dist_thresh)
                    if n > 0:
                        self._interacting_chains.append(x)
                        self._n_pair_contacts.append(n)
        return self._interacting_chains

    @property
    def potentially_contributing_chains(self):
        """ Pairs of chains in :attr:`view` with potential contribution to lDDT

        That are pairs of chains that have at least one interaction within
        :attr:`~dist_thresh` + max(:attr:`~dist_diff_thresholds`)
        """
        if self._potentially_contributing_chains is None:
            self._potentially_contributing_chains = list()
            max_dist_diff_thresh = max(self.dist_diff_thresholds)
            thresh = self.dist_thresh + max_dist_diff_thresh
            for x in itertools.combinations(self.chain_names, 2):
                if self.PotentialInteraction(x[0], x[1],
                                             slack = max_dist_diff_thresh):
                    n = np.count_nonzero(self.PairDist(x[0], x[1]) < thresh)
                    if n > 0:
                        self._potentially_contributing_chains.append(x)

        return self._potentially_contributing_chains

    @property
    def n_pair_contacts(self):
        """ Number of contacts in :attr:`~interacting_chains`

        :type: :class:`list` of :class:`int`
        """
        if self._n_pair_contacts:
            # ugly hack: assumption that computing self.interacting_chains
            # also triggers computation of n_pair_contacts
            int_chains = self.interacting_chains
        return self._n_pair_contacts

    @property
    def n_sc_contacts(self):
        """ Number of contacts for single chains in :attr:`~chain_names`

        :type: :class:`list` of :class:`int`
        """
        if self._n_sc_contacts is None:
            self._n_sc_contacts = list()
            for cname in self.chain_names:
                dist_mat = self.Dist(cname)
                n = np.count_nonzero(dist_mat < self.dist_thresh)
                # dist_mat is symmetric => first remove the diagonal from n
                # as these are distances with itself, i.e. zeroes.
                # Division by two then removes the symmetric component.
                self._n_sc_contacts.append(int((n-dist_mat.shape[0])/2))
        return self._n_sc_contacts

    @property
    def n_contacts(self):
        """ Total number of contacts

        That's the sum of all :attr:`~n_pair_contacts` and
        :attr:`~n_sc_contacts`.

        :type: :class:`int`
        """
        if self._n_contacts is None:
            self._n_contacts = sum(self.n_pair_contacts) + sum(self.n_sc_contacts)
        return self._n_contacts
    
    def GetChain(self, chain_name):
        """ Get chain by name

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """ 
        chain = self.view.FindChain(chain_name)
        if not chain.IsValid():
            raise RuntimeError(f"view has no chain named \"{chain_name}\"")
        return chain

    def GetSequence(self, chain_name):
        """ Get sequence of chain

        Returns sequence of specified chain as raw :class:`str`

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """
        if chain_name not in self._sequence:
            ch = self.GetChain(chain_name)
            s = ''.join([r.one_letter_code for r in ch.residues])
            self._sequence[chain_name] = s
        return self._sequence[chain_name]

    def GetPos(self, chain_name):
        """ Get representative positions of chain

        That's CA positions for peptide residues and C3' for
        nucleotides. Returns positions as a numpy array of shape
        (n_residues, 3).

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """
        if chain_name not in self._pos:
            ch = self.GetChain(chain_name)
            pos = np.zeros((ch.GetResidueCount(), 3))
            for i, r in enumerate(ch.residues):
                pos[i,:] = r.atoms[0].GetPos().data
            self._pos[chain_name] = pos
        return self._pos[chain_name]

    def Dist(self, chain_name):
        """ Get pairwise distance of residues within same chain

        Returns distances as square numpy array of shape (a,a)
        where a is the number of residues in specified chain.
        """
        if chain_name not in self._sc_dist:
            self._sc_dist[chain_name] = distance.cdist(self.GetPos(chain_name),
                                                       self.GetPos(chain_name),
                                                       'euclidean')
        return self._sc_dist[chain_name]

    def PairDist(self, chain_name_one, chain_name_two):
        """ Get pairwise distances between two chains

        Returns distances as numpy array of shape (a, b).
        Where a is the number of residues of the chain that comes BEFORE the
        other in :attr:`~chain_names` 
        """
        key = (min(chain_name_one, chain_name_two),
               max(chain_name_one, chain_name_two))
        if key not in self._pair_dist:
            self._pair_dist[key] = distance.cdist(self.GetPos(key[0]),
                                                  self.GetPos(key[1]),
                                                  'euclidean')
        return self._pair_dist[key]

    def GetMinPos(self, chain_name):
        """ Get min x,y,z cooridnates for given chain

        Based on positions that are extracted with GetPos

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """
        if chain_name not in self._min_pos:
            self._min_pos[chain_name] = self.GetPos(chain_name).min(0)
        return self._min_pos[chain_name]

    def GetMaxPos(self, chain_name):
        """ Get max x,y,z cooridnates for given chain

        Based on positions that are extracted with GetPos

        :param chain_name: Chain in :attr:`~view`
        :type chain_name: :class:`str`
        """
        if chain_name not in self._max_pos:
            self._max_pos[chain_name] = self.GetPos(chain_name).max(0)
        return self._max_pos[chain_name]

    def PotentialInteraction(self, chain_name_one, chain_name_two,
                             slack=0.0):
        """ Returns True if chains potentially interact

        Based on crude collision detection. There is no guarantee
        that they actually interact if True is returned. However,
        if False is returned, they don't interact for sure.

        :param chain_name_one: Chain in :attr:`~view`
        :type chain_name_one: class:`str`
        :param chain_name_two: Chain in :attr:`~view`
        :type chain_name_two: class:`str`
        :param slack: Add slack to interaction distance threshold
        :type slack: :class:`float`
        """
        min_one = self.GetMinPos(chain_name_one)
        max_one = self.GetMaxPos(chain_name_one)
        min_two = self.GetMinPos(chain_name_two)
        max_two = self.GetMaxPos(chain_name_two)
        if np.max(min_one - max_two) > (self.dist_thresh + slack):
            return False
        if np.max(min_two - max_one) > (self.dist_thresh + slack):
            return False
        return True


class BBlDDTScorer:
    """ Helper object to compute Backbone only lDDT score

    Tightly integrated into the mechanisms from the chain_mapping module.
    The prefered way to derive an object of type :class:`BBlDDTScorer` is
    through the static constructor: :func:`~FromMappingResult`.

    lDDT computation in :func:`BBlDDTScorer.Score` implements caching.
    Repeated computations with alternative chain mappings thus become faster.

    :param target: Structure designated as "target". Can be fetched from
                   :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type target: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param chem_groups: Groups of chemically equivalent chains in *target*.
                        Can be fetched from
                        :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type chem_groups: :class:`list` of :class:`list` of :class:`str`
    :param model: Structure designated as "model". Can be fetched from
                  :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param alns: Each alignment is accessible with ``alns[(t_chain,m_chain)]``.
                 First sequence is the sequence of the respective chain in
                 :attr:`~qsent1`, second sequence the one from :attr:`~qsent2`.
                 Can be fetched from
                 :class:`ost.mol.alg.chain_mapping.MappingResult`
    :type alns: :class:`dict` with key: :class:`tuple` of :class:`str`, value:
                :class:`ost.seq.AlignmentHandle`
    :param dist_thresh: Max distance of a pairwise interaction in target
                        to be considered as contact in lDDT
    :type dist_thresh: :class:`float`
    :param dist_diff_thresholds: Distance difference thresholds for
                                 lDDT computations
    :type dist_diff_thresholds: :class:`list` of :class:`float`
    """
    def __init__(self, target, chem_groups, model, alns, dist_thresh = 15.0,
                 dist_diff_thresholds = [0.5, 1.0, 2.0, 4.0]):

        self._trg = BBlDDTEntity(target, dist_thresh = dist_thresh,
                                 dist_diff_thresholds=dist_diff_thresholds)

        # ensure that target chain names match the ones in chem_groups
        chem_group_ch_names = list(itertools.chain.from_iterable(chem_groups))
        if self._trg.chain_names != sorted(chem_group_ch_names):
            raise RuntimeError(f"Expect exact same chain names in chem_groups "
                               f"and in target (which is processed to only "
                               f"contain peptides/nucleotides). target: "
                               f"{self._trg.chain_names}, chem_groups: "
                               f"{chem_group_ch_names}")

        self._chem_groups = chem_groups
        self._mdl = BBlDDTEntity(model, dist_thresh = dist_thresh,
                                 dist_diff_thresholds=dist_diff_thresholds)
        self._alns = alns
        self._dist_diff_thresholds = dist_diff_thresholds
        self._dist_thresh = dist_thresh

        # cache for mapped interface scores
        # key: tuple of tuple ((trg_ch1, trg_ch2),
        #                     ((mdl_ch1, mdl_ch2))
        # where the first tuple is lexicographically sorted
        # the second tuple is mapping dependent
        # value: numpy array of len(dist_thresholds) representing the
        # respective numbers of fulfilled contacts
        self._pairwise_cache = dict()

        # cache for mapped single chain scores
        # key: tuple (trg_ch, mdl_ch)
        # value: numpy array of len(dist_thresholds) representing the
        # respective numbers of fulfilled contacts
        self._sc_cache = dict()

    @staticmethod
    def FromMappingResult(mapping_result, dist_thresh = 15.0,
                          dist_diff_thresholds = [0.5, 1.0, 2.0, 4.0]):
        """ The preferred way to get a :clas:`BBlDDTScorer`

        Static constructor that derives an object of type :class:`QSScorer`
        using a :class:`ost.mol.alg.chain_mapping.MappingResult`

        :param mapping_result: Data source
        :type mapping_result: :class:`ost.mol.alg.chain_mapping.MappingResult`
        :param dist_thresh: The lDDT distance threshold
        :type dist_thresh: :class:`float`
        :param dist_diff_thresholds: The lDDT distance difference thresholds
        :type dist_diff_thresholds: :class:`list` of :class:`float`       
        """
        scorer = BBlDDTScorer(mapping_result.target, mapping_result.chem_groups,
                              mapping_result.model, alns = mapping_result.alns,
                              dist_thresh = dist_thresh,
                              dist_diff_thresholds = dist_diff_thresholds)
        return scorer

    @property
    def trg(self):
        """ The :class:`BBlDDTEntity` representing target

        :type: :class:`BBlDDTEntity`
        """
        return self._trg

    @property
    def mdl(self):
        """ The :class:`BBlDDTEntity` representing model

        :type: :class:`BBlDDTEntity`
        """
        return self._mdl

    @property
    def alns(self):
        """ Alignments between chains in :attr:`~trg` and :attr:`~mdl`

        Provided at object construction. Each alignment is accessible with
        ``alns[(t_chain,m_chain)]``. First sequence is the sequence of the
        respective chain in :attr:`~trg`, second sequence the one from
        :attr:`~mdl`.

        :type: :class:`dict` with key: :class:`tuple` of :class:`str`, value:
               :class:`ost.seq.AlignmentHandle`
        """
        return self._alns

    @property
    def chem_groups(self):
        """ Groups of chemically equivalent chains in :attr:`~trg`

        Provided at object construction

        :type: :class:`list` of :class:`list` of :class:`str`
        """
        return self._chem_groups

    def Score(self, mapping, check=True):
        """ Computes Backbone lDDT given chain mapping

        Again, the preferred way is to get *mapping* is from an object
        of type :class:`ost.mol.alg.chain_mapping.MappingResult`.

        :param mapping: see 
                        :attr:`ost.mol.alg.chain_mapping.MappingResult.mapping`
        :type mapping: :class:`list` of :class:`list` of :class:`str`
        :param check: Perform input checks, can be disabled for speed purposes
                      if you know what you're doing.
        :type check: :class:`bool`
        :returns: The score
        """
        if check:
            # ensure that dimensionality of mapping matches self.chem_groups
            if len(self.chem_groups) != len(mapping):
                raise RuntimeError("Dimensions of self.chem_groups and mapping "
                                   "must match")
            for a,b in zip(self.chem_groups, mapping):
                if len(a) != len(b):
                    raise RuntimeError("Dimensions of self.chem_groups and "
                                       "mapping must match")
            # ensure that chain names in mapping are all present in qsent2
            for name in itertools.chain.from_iterable(mapping):
                if name is not None and name not in self.mdl.chain_names:
                    raise RuntimeError(f"Each chain in mapping must be present "
                                       f"in self.mdl. No match for "
                                       f"\"{name}\"")

        flat_mapping = dict()
        for a, b in zip(self.chem_groups, mapping):
            flat_mapping.update({x: y for x, y in zip(a, b) if y is not None})

        return self.FromFlatMapping(flat_mapping)

    def FromFlatMapping(self, flat_mapping):
        """ Same as :func:`Score` but with flat mapping

        :param flat_mapping: Dictionary with target chain names as keys and
                             the mapped model chain names as value
        :type flat_mapping: :class:`dict` with :class:`str` as key and value
        :returns: :class:`float` representing lDDT
        """
        n_conserved = np.zeros(len(self._dist_diff_thresholds), dtype=np.int32)

        # process single chains
        for cname in self.trg.chain_names:
            if cname in flat_mapping:
                n_conserved += self._NSCConserved(cname, flat_mapping[cname])

        # process interfaces
        for interface in self.trg.interacting_chains:
            if interface[0] in flat_mapping and interface[1] in flat_mapping:
                mdl_interface = (flat_mapping[interface[0]],
                                 flat_mapping[interface[1]])
                n_conserved += self._NPairConserved(interface, mdl_interface)

        return np.mean(n_conserved / self.trg.n_contacts)

    def _NSCConserved(self, trg_ch, mdl_ch):
        if (trg_ch, mdl_ch) in self._sc_cache:
            return self._sc_cache[(trg_ch, mdl_ch)]
        trg_dist = self.trg.Dist(trg_ch)
        mdl_dist = self.mdl.Dist(mdl_ch)
        trg_indices, mdl_indices = self._IndexMapping(trg_ch, mdl_ch)
        trg_dist = trg_dist[np.ix_(trg_indices, trg_indices)]
        mdl_dist = mdl_dist[np.ix_(mdl_indices, mdl_indices)]
        # mask to select relevant distances (dist in trg < dist_thresh)
        # np.triu zeroes the values below the diagonal
        mask = np.triu(trg_dist < self._dist_thresh)
        n_diag = trg_dist.shape[0]
        trg_dist = trg_dist[mask]
        mdl_dist = mdl_dist[mask]
        dist_diffs = np.absolute(trg_dist - mdl_dist)
        n_conserved = np.zeros(len(self._dist_diff_thresholds), dtype=np.int32)
        for thresh_idx, thresh in enumerate(self._dist_diff_thresholds):
            N = (dist_diffs < thresh).sum()
            # still need to consider the 0.0 dist diffs on the diagonal
            n_conserved[thresh_idx] = int((N - n_diag))
        self._sc_cache[(trg_ch, mdl_ch)] = n_conserved
        return n_conserved

    def _NPairConserved(self, trg_int, mdl_int):
        key_one = (trg_int, mdl_int)
        if key_one in self._pairwise_cache:
            return self._pairwise_cache[key_one]
        key_two = ((trg_int[1], trg_int[0]), (mdl_int[1], mdl_int[0]))
        if key_two in self._pairwise_cache:
            return self._pairwise_cache[key_two]
        trg_dist = self.trg.PairDist(trg_int[0], trg_int[1])
        mdl_dist = self.mdl.PairDist(mdl_int[0], mdl_int[1])
        if trg_int[0] > trg_int[1]:
            trg_dist = trg_dist.transpose()
        if mdl_int[0] > mdl_int[1]:
            mdl_dist = mdl_dist.transpose()
        trg_indices_1, mdl_indices_1 = self._IndexMapping(trg_int[0], mdl_int[0])
        trg_indices_2, mdl_indices_2 = self._IndexMapping(trg_int[1], mdl_int[1])
        trg_dist = trg_dist[np.ix_(trg_indices_1, trg_indices_2)]
        mdl_dist = mdl_dist[np.ix_(mdl_indices_1, mdl_indices_2)]
        # reduce to relevant distances (dist in trg < dist_thresh)
        mask = trg_dist < self._dist_thresh
        trg_dist = trg_dist[mask]
        mdl_dist = mdl_dist[mask]
        dist_diffs = np.absolute(trg_dist - mdl_dist)
        n_conserved = np.zeros(len(self._dist_diff_thresholds), dtype=np.int32)
        for thresh_idx, thresh in enumerate(self._dist_diff_thresholds):
            n_conserved[thresh_idx] = (dist_diffs < thresh).sum()
        self._pairwise_cache[key_one] = n_conserved
        return n_conserved

    def _IndexMapping(self, ch1, ch2):
        """ Fetches aln and returns indices of aligned residues

        returns 2 numpy arrays containing the indices of residues in
        ch1 and ch2 which are aligned
        """
        mapped_indices_1 = list()
        mapped_indices_2 = list()
        idx_1 = 0
        idx_2 = 0
        for col in self.alns[(ch1, ch2)]:
            if col[0] != '-' and col[1] != '-':
                mapped_indices_1.append(idx_1)
                mapped_indices_2.append(idx_2)
            if col[0] != '-':
                idx_1 +=1
            if col[1] != '-':
                idx_2 +=1
        return (np.array(mapped_indices_1), np.array(mapped_indices_2))
