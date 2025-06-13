import ost
from ost import io
from ost import conop
from ost import mol
from ost import seq


def CleanHydrogens(ent, clib):
    """ Scoring helper - Returns copy of *ent* without hydrogens

    Non-standard hydrogen naming can cause trouble in residue property
    assignment which is done by the :class:`ost.conop.RuleBasedProcessor` when
    loading. In fact, residue property assignment is not done for every residue
    that has unknown atoms according to the chemical component dictionary. This
    function therefore re-processes the entity after removing hydrogens.

    :param ent: Entity to clean
    :type ent: :class:`ost.mol.EntityHandle`/:class:`ost.mol.EntityView`
    :param clib: Compound library to perform re-processing after hydrogen
                 removal.
    :type clib: :class:`ost.conop.CompoundLib`
    :returns: Cleaned and re-processed ent
    """
    cleaned_ent = mol.CreateEntityFromView(ent.Select(
        "ele != H and ele != D"), include_exlusive_atoms=False)
    # process again to set missing residue properties due to non standard
    # hydrogens
    processor = conop.RuleBasedProcessor(clib)
    processor.Process(cleaned_ent)
    return cleaned_ent


def MMCIFPrep(mmcif_path, biounit=None, extract_nonpoly=False,
              fault_tolerant=False, extract_seqres_mapping=False):
    """ Scoring helper - Prepares input from mmCIF

    Only performs gentle cleanup of hydrogen atoms. Further cleanup is delegated
    to scoring classes.

    Depending on input flags, the following outputs can be retrieved:

    * poly_ent (:class:`ost.mol.EntityHandle`): An OpenStructure entity with only
      polymer chains. This is based on _entity.type extracted from *mmcif_path*.
      If _entity.type is not defined for every chain, a warning is logged
      and the returned poly_ent is a selection for peptide and nucleotide
      residues as defined in the chemical component dictionary.
    * non_poly_entities (:class:`list` of :class:`ost.mol.EntityHandle`):
      OpenStructure entities representing all non-polymer (ligand) entities.
      This is based on _entity.type extracted from *mmcif_path*. If _entity.type
      is not defined for every chain, an error is raised.
    * seqres (:class:`ost.seq.SequenceList`): Seqres sequences with entity id
      as sequence names and the respective canonical seqres as sequence. Set to
      None and triggers a warning if information is missing in *mmcif_path*.
    * trg_seqres_mapping (:class:`dict`): Dictionary with chain names in
      poly_ent as keys and the respective entity ids as values.
      Set to None and triggers a warning if information is missing in
      *mmcif_path*.

    :param mmcif_path: Path to mmCIF file that contains polymer and optionally
                       non-polymer entities
    :type mmcif_path: :class:`str`
    :param biounit: If given, construct specified biounit from mmCIF AU
    :type biounit: :class:`str`
    :param extract_nonpoly: Controls return value
    :type extract_nonpoly: :class:`bool`
    :param fault_tolerant: Passed as parameter to :func:`ost.io.LoadMMCIF`
    :type fault_tolerant: :class:`bool`
    :param extract_seqres_mapping: Controls return value
    :type extract_seqres_mapping: :class:`bool`
    :returns: poly_ent if *extract_nonpoly*/*extract_seqres_mapping* are False.
              (poly_ent, non_poly_entities) if *extract_nonpoly* is True.
              (poly_ent, seqres, trg_seqres_mapping) if *extract_seqres_mapping*
              is True.
              (poly_ent, non_poly_entities, seqres, trg_seqres_mapping) if both
              flags are True.
    """
    clib = conop.GetDefaultLib()
    if not clib:
        ost.LogError("A compound library is required. "
                     "Please refer to the OpenStructure website: "
                     "https://openstructure.org/docs/conop/compoundlib/.")
        raise RuntimeError("No compound library found")

    # return variables that will be defined depending on input flags
    poly_ent = None
    non_poly_entities = None
    seqres = None
    trg_seqres_mapping = None

    # increase loglevel, as we would pollute the info log with weird stuff
    ost.PushVerbosityLevel(ost.LogLevel.Error)
    mmcif_entity, mmcif_seqres, mmcif_info = io.LoadMMCIF(mmcif_path, seqres=True, info=True,
                                                          fault_tolerant=fault_tolerant)
    # restore old loglevel and return
    ost.PopVerbosityLevel()

    mmcif_entity = CleanHydrogens(mmcif_entity, clib)

    # construct biounit if necessary
    if biounit is not None:
        biounit_found = False
        for bu in mmcif_info.biounits:
            if bu.id == biounit:
                mmcif_entity = mol.alg.CreateBU(mmcif_entity, bu)
                biounit_found = True
                break
        if not biounit_found:
            raise RuntimeError(f"Specified biounit '{biounit}' not in "
                               f"{mmcif_path}")

    # check if we have entity types defined for each chain
    missing_entity_types = list()
    for ch in mmcif_entity.chains:
        cname = None
        if biounit is not None:
            # if a biounit is constructed, you get chain names like: 1.YOLO
            # we cannot simply split by '.' since '.' is an allowed character
            # in chain names. => split by first occurence
            dot_index = ch.name.find('.')
            if dot_index == -1:
                cname = ch.name
            else:
                cname = ch.name[dot_index+1:]
        else:
            cname = ch.name
        try:
            entity_id = mmcif_info.GetMMCifEntityIdTr(cname)
            # the following raises if there is no desc for entity_id
            entity_desc = mmcif_info.GetEntityDesc(entity_id)
        except:
            missing_entity_types.append(cname)

    if len(missing_entity_types) > 0:
        msg = f"mmCIF file does not define _entity.type for chains "
        msg += f"{missing_entity_types}. "
        
        if not fault_tolerant:
            msg += f"Use fault tolerant mode to ignore this error and "
            msg += f"fallback to select polymers based on peptide/nucleotide "
            msg += f"residues as defined in chemical component dictionary."
            raise RuntimeError(msg)

        msg += f"Fallback to select polymers based on peptide/nucleotide "
        msg += f"residues as defined in chemical component dictionary (fault "
        msg += f"tolerant mode)."
        ost.LogWarning(msg)

        poly_sel = mmcif_entity.Select("peptide=true or nucleotide=true")
        poly_ent = mol.CreateEntityFromView(poly_sel, True)
    else:
        polymer_entity_ids = mmcif_info.GetEntityIdsOfType("polymer")
        for ch in mmcif_entity.chains:
            cname = None
            if biounit is not None:
                # if a biounit is constructed, you get chain names like: 1.YOLO
                # we cannot simply split by '.' since '.' is an allowed character
                # in chain names. => split by first occurence
                dot_index = ch.name.find('.')
                if dot_index == -1:
                    cname = ch.name
                else:
                    cname = ch.name[dot_index+1:]
            else:
                cname = ch.name
            if mmcif_info.GetMMCifEntityIdTr(cname) in polymer_entity_ids:
                ch.SetIntProp("poly", 1)
        poly_sel = mmcif_entity.Select("gcpoly:0=1")
        poly_ent = mol.CreateEntityFromView(poly_sel, True)

    if extract_nonpoly:
        if len(missing_entity_types) > 0:
            msg = f"mmCIF file does not contain _entity.type for the following "
            msg += f"chain(s): {missing_entity_types}. Extracting non-polymers "
            msg += f"from mmCIF files requires _entity.type to be set for all "
            msg += f"chains."
            raise RuntimeError(msg)

        non_polymer_entity_ids = mmcif_info.GetEntityIdsOfType("non-polymer")
        nonpoly_id = 1
        for ch in mmcif_entity.chains:
            cname = None
            if biounit is not None:
                # if a biounit is constructed, you get chain names like: 1.YOLO
                # we cannot simply split by '.' since '.' is an allowed character
                # in chain names. => split by first occurence
                dot_index = ch.name.find('.')
                if dot_index == -1:
                    cname = ch.name
                else:
                    cname = ch.name[dot_index+1:]
            else:
                cname = ch.name
            if mmcif_info.GetMMCifEntityIdTr(cname) in non_polymer_entity_ids:
                ch.SetIntProp("nonpolyid", nonpoly_id)
                nonpoly_id += 1

        non_poly_entities = list()
        for i in range(1,nonpoly_id):
            view = mmcif_entity.Select(f"gcnonpolyid:0={i}")
            if view.GetResidueCount() != 1:
                raise RuntimeError(f"Expected non-polymer entities in "
                                   f"{mmcif_path} to contain exactly 1 "
                                   f"residue. Got {view.GetResidueCount()} "
                                   f"in chain {view.chains[0].name}")
            compound = clib.FindCompound(view.residues[0].name)
            if compound is None:
                error_msg = f"\"{view.residues[0].name}\" is not available in " \
                            f"the compound library."
                if fault_tolerant:
                    error_msg += f"A distance-based heuristic was used to " \
                                 f"connect the ligand atoms (fault tolerant " \
                                 f"mode)."
                    ost.LogWarning(error_msg)
                else:
                    error_msg += f"Use fault tolerant mode to ignore this " \
                                 f"error and use a distance based heuristic " \
                                 f"to connect the ligand atoms."
                    raise RuntimeError(error_msg)

            non_poly_entities.append(mol.CreateEntityFromView(view, True))

    if extract_seqres_mapping:
        # mmcif seqres is a list of sequences that relates to
        # chain names in the assymetric unit. What we want is a list
        # of sequences that relate to the underlying entities.
        seqres = seq.CreateSequenceList()
        seqres_processed = set()

        for s in mmcif_seqres:
            entity_id = mmcif_info.GetMMCifEntityIdTr(s.GetName())
            if entity_id not in seqres_processed:
                seqres_processed.add(entity_id)
                seqres.AddSequence(seq.CreateSequence(entity_id, s.GetGaplessString()))

        # check if we have SEQRES defined for each polymer chain
        missing_seqres = list()
        for ch in poly_ent.chains:
            cname = None
            if biounit is not None:
                # if a biounit is constructed, you get chain names like: 1.YOLO
                # we cannot simply split by '.' since '.' is an allowed character
                # in chain names. => split by first occurence
                dot_index = ch.name.find('.')
                if dot_index == -1:
                    cname = ch.name
                else:
                    cname = ch.name[dot_index+1:]
            else:
                cname = ch.name
            
            entity_id = mmcif_info.GetMMCifEntityIdTr(cname)
            if entity_id not in seqres_processed:
                missing_seqres.append(cname)

        if len(missing_seqres) > 0:
            msg = f"Extracting chem grouping from mmCIF file requires all "
            msg += f"SEQRES information set. SEQRES is missing for polymer "
            msg += f"chain(s) {missing_seqres}. "
            
            if not fault_tolerant:
                msg += f"Use fault tolerant mode to ignore this error and "
                msg += f"fallback sequence identity-based chem grouping."
                raise RuntimeError(msg)
            
            msg += f"Chem grouping will be based on sequence identity (fault "
            msg += f"tolerant mode)."
            ost.LogWarning(msg)

            seqres = None
            trg_seqres_mapping = None
        else:
            trg_seqres_mapping = dict()
            if biounit is None:
                cnames = [ch.name for ch in poly_ent.chains]
                for cname in cnames:
                    trg_seqres_mapping[cname] = mmcif_info.GetMMCifEntityIdTr(cname)
            else:
                bu_cnames = [ch.name for ch in poly_ent.chains]
                au_cnames = list()
                for bu_cname in bu_cnames:
                    dot_idx = bu_cname.index(".")
                    au_cnames.append(bu_cname[dot_idx + 1 :])
                for au_cname, bu_cname in zip(au_cnames, bu_cnames):
                    trg_seqres_mapping[bu_cname] = mmcif_info.GetMMCifEntityIdTr(au_cname)

    if extract_nonpoly and extract_seqres_mapping:
        return (poly_ent, non_poly_entities, seqres, trg_seqres_mapping)
    elif extract_nonpoly:
        return (poly_ent, non_poly_entities)
    elif extract_seqres_mapping:
        return (poly_ent, seqres, trg_seqres_mapping)
    else:
        return poly_ent


def PDBPrep(pdb_path, fault_tolerant=False):
    """ Scoring helper - Prepares scoring input from PDB

    Only performs gentle cleanup of hydrogen atoms. Further cleanup is delegated
    to scoring classes. There is no logic to extract ligands from PDB
    files. Ligands must be provided separately as SDF files in these cases.

    :param pdb_path: Path to PDB file that contains polymer entities
    :type pdb_path: :class:`str`
    :param fault_tolerant: Passed as parameter to :func:`ost.io.LoadPDB`
    :type fault_tolerant: :class:`bool`
    :returns: :class:`EntityHandle` from loaded file.
    """
    clib = conop.GetDefaultLib()
    if not clib:
        ost.LogError("A compound library is required. "
                     "Please refer to the OpenStructure website: "
                     "https://openstructure.org/docs/conop/compoundlib/.")
        raise RuntimeError("No compound library found")

    # increase loglevel, as we would pollute the info log with weird stuff
    ost.PushVerbosityLevel(ost.LogLevel.Error)
    pdb_entity = io.LoadPDB(pdb_path, fault_tolerant=fault_tolerant)
    # restore old loglevel and return
    ost.PopVerbosityLevel()
    pdb_entity = CleanHydrogens(pdb_entity, clib)

    return pdb_entity

__all__ = ('CleanHydrogens', 'MMCIFPrep', 'PDBPrep')
