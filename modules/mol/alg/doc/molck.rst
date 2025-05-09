.. _molck:

Molecular Checker (Molck)
--------------------------------------------------------------------------------

The Molecular Checker (Molck) is a tool for cleaning up molecular structures 
and making them conform to the :doc:`compound library  <../../conop/compoundlib>`.

Molck removes any residues and atoms that are not defined in the compound
library. This means that if the structure contains residues or atoms that 
are not part of the compound library, they will be removed during the cleaning
process.

.. caution::
  Do not use Molck if you need to preserve residues or atoms that are not
  defined in the compound library. For example, if your structure contains
  ligands or other custom molecules that are not in the compound library,
  using Molck would not preserve these components.

Programmatic usage
##################

Molecular Checker (Molck) could be called directly from the code using Molck
function:

.. code-block:: python

  #! /bin/env python

  """Run Molck with Python API.


  This is an exemplary procedure on how to run Molck using Python API which is
  equivalent to the command line:

  molck <PDB PATH> --rm=hyd,oxt,nonstd,unk \
                   --fix-ele --out=<OUTPUT PATH> \
                   --complib=<PATH TO compounds.chemlib>
  """

  from ost.io import LoadPDB, SavePDB
  from ost.mol.alg import MolckSettings, Molck
                         
  from ost.conop import CompoundLib


  pdbid = "<PDB PATH>"
  lib = CompoundLib.Load("<PATH TO compounds.chemlib>")

  # Using Molck function
  ent = LoadPDB(pdbid)
  ms = MolckSettings(rm_unk_atoms=True,
                     rm_non_std=True,
                     rm_hyd_atoms=True,
                     rm_oxt_atoms=True,
                     rm_zero_occ_atoms=False,
                     colored=False,
                     map_nonstd_res=False,
                     assign_elem=True)
  Molck(ent, lib, ms)
  SavePDB(ent, "<OUTPUT PATH>")

It can also be split into subsequent commands for greater controll:

.. code-block:: python

  #! /bin/env python

  """Run Molck with Python API.


  This is an exemplary procedure on how to run Molck using Python API which is
  equivalent to the command line:

  molck <PDB PATH> --rm=hyd,oxt,nonstd,unk \
                   --fix-ele --out=<OUTPUT PATH> \
                   --complib=<PATH TO compounds.chemlib>
  """

  from ost.io import LoadPDB, SavePDB
  from ost.mol.alg import (RemoveAtoms, MapNonStandardResidues,
                           CleanUpElementColumn)
  from ost.conop import CompoundLib


  pdbid = "<PDB PATH>"
  lib = CompoundLib.Load("<PATH TO compounds.chemlib>")
  map_nonstd = False

  # Using function chain
  ent = LoadPDB(pdbid)
  if map_nonstd:
      MapNonStandardResidues(lib=lib, ent=ent)

  RemoveAtoms(lib=lib,
              ent=ent,
              rm_unk_atoms=True,
              rm_non_std=True,
              rm_hyd_atoms=True,
              rm_oxt_atoms=True,
              rm_zero_occ_atoms=False,
              colored=False)

  CleanUpElementColumn(lib=lib, ent=ent)
  SavePDB(ent, "<OUTPUT PATH>")

API
###

.. class:: MolckSettings(rm_unk_atoms=True, rm_non_std=False, \
                         rm_hyd_atoms=True, rm_oxt_atoms=False, \
                         rm_zero_occ_atoms=False, colored=False, \
                         map_nonstd_res=True, assign_elem=True)

  Stores settings used for Molecular Checker.

  :param rm_unk_atoms: Sets :attr:`rm_unk_atoms`.
  :param rm_non_std: Sets :attr:`rm_non_std`.
  :param rm_hyd_atoms: Sets :attr:`rm_hyd_atoms`.
  :param rm_oxt_atoms: Sets :attr:`rm_oxt_atoms`.
  :param rm_zero_occ_atoms: Sets :attr:`rm_zero_occ_atoms`.
  :param colored: Sets :attr:`colored`.
  :param map_nonstd_res: Sets :attr:`map_nonstd_res`.
  :param assign_elem: Sets :attr:`assign_elem`.

  .. attribute:: rm_unk_atoms

    .. tip::

      This flag should **always** be set to True. Other flags will behave
      unexpectedly otherwise.

    Remove unknown atoms. That is 1) any atom from residues that are not
    present in the compound library (provided at Molck call) and 2) any atom
    with a name that is not present in the respective entries of the compound
    library.
    
    :type: :class:`bool`

  .. attribute:: rm_non_std

    Remove all residues not one of the 20 standard amino acids.
    This removes all other residues including unknown residues, ligands,
    saccharides and nucleotides (including the 4 standard nucleotides).
    
    :type: :class:`bool`

  .. attribute:: rm_hyd_atoms

    Remove hydrogen atoms. That's all atoms with element specified as H or D
    in the respective entries of the compound library (provided at Molck call).
    Unknown atoms (see :attr:`rm_unk_atoms`) are not removed by this flag. If you
    really want to get rid of every hydrogen, you need to combine it with
    :attr:`rm_unk_atoms`.
    
    :type: :class:`bool`

  .. attribute:: rm_oxt_atoms

    Remove all atoms with name "OXT". That's typically terminal oxygens in protein
    chains, but this might remove arbitrary atoms in other molecules. You should
    only use this flag in combination with :attr:`rm_non_std`.
    
    :type: :class:`bool`

  .. attribute:: rm_zero_occ_atoms

    Remove atoms with zero occupancy.
    
    :type: :class:`bool`

  .. attribute:: colored

    Whether output should be colored.
    
    :type: :class:`bool`

  .. attribute:: map_nonstd_res

    Maps modified residues back to the parent amino acid, for example
    MSE -> MET, SEP -> SER.
    
    :type: :class:`bool`

  .. attribute:: assign_elem

    Assigns elements as defined in the respective entries of the compound 
    library (provided at Molck call). For unknown atoms (see definition in
    :attr:`rm_unk_atoms`), the element is set to an empty string.
    To avoid empty strings as elements, this property should only be applied
    in combination with :attr:`rm_unk_atoms`.
    
    :type: :class:`bool`

  .. method:: ToString()

    :return: String representation of the MolckSettings.
    :rtype:  :class:`str`

.. warning::

  The API here is set such that the functions modify the passed structure *ent*
  in-place. If this is not ok, please work on a copy of the structure.

.. function:: Molck(ent, lib, settings, [prune=True])

  Runs Molck on provided entity. Reprocesses *ent* with
  :class:`ost.conop.HeuristicProcessor` and given *lib* once done.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`
  :param settings: Molck settings
  :type settings: :class:`MolckSettings`
  :param prune: Whether to remove residues/chains that don't contain atoms 
                anymore after Molck cleanup
  :type prune: :class:`bool` 


.. function:: MapNonStandardResidues(ent, lib, reprocess=True)

  Maps modified residues back to the parent amino acid, for example MSE -> MET.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`
  :param reprocess: The function generates a deep copy of *ent*. Highly
                    recommended to enable *reprocess* that runs
                    :class:`ost.conop.HeuristicProcessor` with given *lib*.
                    If set to False, you'll have no connectivity etc. after
                    calling this function.

.. function:: RemoveAtoms(ent, lib, rm_unk_atoms=True, rm_non_std=False, \
                          rm_hyd_atoms=True, rm_oxt_atoms=False, \
                          rm_zero_occ_atoms=False, colored=False,
                          reprocess=True)

  Removes atoms and residues according to some criteria.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`
  :param rm_unk_atoms: See :attr:`MolckSettings.rm_unk_atoms`
  :param rm_non_std: See :attr:`MolckSettings.rm_non_std`
  :param rm_hyd_atoms: See :attr:`MolckSettings.rm_hyd_atoms`
  :param rm_oxt_atoms: See :attr:`MolckSettings.rm_oxt_atoms`
  :param rm_zero_occ_atoms: See :attr:`MolckSettings.rm_zero_occ_atoms`
  :param colored: See :attr:`MolckSettings.colored`
  :param reprocess: Removing atoms may impact certain annotations on the
                    structure (chem class etc.) which are set by 
                    :class:`ost.conop.Processor`. If set to True,
                    a :class:`ost.conop.HeuristicProcessor` with given
                    *lib* reprocesses *ent*.

.. function:: CleanUpElementColumn(ent, lib)

  Assigns elements as defined in the respective entries of the compound library
  as described in :attr:`MolckSettings.assign_elem`. This should only be called
  after :func:`RemoveAtoms` with :attr:`rm_unk_atoms` set to True.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`