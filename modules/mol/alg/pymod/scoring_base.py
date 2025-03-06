import ost
from ost import io
from ost import conop
from ost import mol


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
              fault_tolerant=False, allow_heuristic_conn=False):
    """ Scoring helper - Prepares input from mmCIF

    Only performs gentle cleanup of hydrogen atoms. Further cleanup is delegated
    to scoring classes.

    :param mmcif_path: Path to mmCIF file that contains polymer and optionally
                       non-polymer entities
    :type mmcif_path: :class:`str`
    :param biounit: If given, construct specified biounit from mmCIF AU
    :type biounit: :class:`str`
    :param extract_nonpoly: Additionally returns a list of
                            :class:`ost.mol.EntityHandle`
                            objects representing all non-polymer (ligand)
                            entities.
    :type extract_nonpoly: :class:`bool`
    :param fault_tolerant: Passed as parameter to :func:`ost.io.LoadMMCIF`
    :type fault_tolerant: :class:`bool`
    :param allow_heuristic_conn: Only relevant if extract_nonpoly is True.
                                 The chemical component dictionary is relevant
                                 for connectivity information. By default, we
                                 enforce the presence of each non-polymer in
                                 the dictionary to ensure correct connectity.
                                 If you enable this flag, you allow the use
                                 of a distance based heuristic as fallback.
                                 With all its consequences in ligand matching.
    :type allow_heuristic_conn: :class:`bool`
    :returns: :class:`ost.mol.EntityHandle` which only contains polymer
              entities representing the receptor structure. If *extract_nonpoly*
              is True, a tuple is returned which additionally contains a
              :class:`list` of :class:`ost.mol.EntityHandle`, where each
              :class:`ost.mol.EntityHandle` represents a non-polymer (ligand).
    """
    clib = conop.GetDefaultLib()
    if not clib:
        ost.LogError("A compound library is required. "
                     "Please refer to the OpenStructure website: "
                     "https://openstructure.org/docs/conop/compoundlib/.")
        raise RuntimeError("No compound library found")

    mmcif_entity, mmcif_info = io.LoadMMCIF(mmcif_path, info=True,
                                            fault_tolerant=fault_tolerant)
    mmcif_entity = CleanHydrogens(mmcif_entity, clib)

    # get AU chain names representing polymer entities
    polymer_entity_ids = mmcif_info.GetEntityIdsOfType("polymer")
    polymer_chain_names = list()
    for ch in mmcif_entity.chains:
        if mmcif_info.GetMMCifEntityIdTr(ch.name) in polymer_entity_ids:
            polymer_chain_names.append(ch.name)

    # get AU chain names representing non-polymer entities
    non_polymer_entity_ids = mmcif_info.GetEntityIdsOfType("non-polymer")
    non_polymer_chain_names = list()
    for ch in mmcif_entity.chains:
        if mmcif_info.GetMMCifEntityIdTr(ch.name) in non_polymer_entity_ids:
            non_polymer_chain_names.append(ch.name)

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

    # assign generic properties for selection later on
    non_poly_id = 0
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
        
        if cname in polymer_chain_names:
            ch.SetIntProp("poly", 1)
        if cname in non_polymer_chain_names:
            ch.SetIntProp("nonpolyid", non_poly_id)
            non_poly_id += 1

    poly_sel = mmcif_entity.Select("gcpoly:0=1")
    poly_ent = mol.CreateEntityFromView(poly_sel, True)

    if extract_nonpoly == False:
        return poly_ent

    non_poly_sel = mmcif_entity.Select("gcnonpoly:0=1")
    non_poly_entities = list()
    for i in range(non_poly_id):
        view = mmcif_entity.Select(f"gcnonpolyid:{non_poly_id}={i}")
        if view.GetResidueCount() != 1:
            raise RuntimeError(f"Expect non-polymer entities in "
                               f"{mmcif_path} to contain exactly 1 "
                               f"residue. Got {ch.GetResidueCount()} "
                               f"in chain {ch.name}")
        if not allow_heuristic_conn:
            compound = clib.FindCompound(view.residues[0].name)
            if compound is None:
                raise RuntimeError(f"Can only extract non-polymer entities if "
                                   f"respective residues are available in PDB "
                                   f"component dictionary. Can't find "
                                   f"\"{view.residues[0].name}\"")

        non_poly_entities.append(mol.CreateEntityFromView(view, True))

    return (poly_ent, non_poly_entities)


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

    pdb_entity = io.LoadPDB(pdb_path, fault_tolerant=fault_tolerant)
    pdb_entity = CleanHydrogens(pdb_entity, clib)

    return pdb_entity

__all__ = ('CleanHydrogens', 'MMCIFPrep', 'PDBPrep')
