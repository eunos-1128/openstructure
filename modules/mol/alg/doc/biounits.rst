Biounits
================================================================================

Biological assemblies, i.e. biounits, are an integral part of mmCIF files and
their construction is fully defined in :class:`ost.io.MMCifInfoBioUnit`.
:func:`ost.io.MMCifInfoBioUnit.PDBize` provides one possibility to construct
such biounits with compatibility with the PDB format in mind. That is single
character chain names, dumping all ligands in one chain etc. Here we provide a
more mmCIF-style way of constructing biounits. This can either be done starting
from a :class:`ost.io.MMCifInfoBioUnit` or the derived
:class:`ost.mol.alg.BUInfo`. The latter is a minimalistic representation of
:class:`ost.io.MMCifInfoBioUnit` and can be serialized to a byte string.

.. class:: BUInfo(mmcif_buinfo):

  Preprocesses data from :class:`ost.io.MMCifInfoBioUnit` that are required
  to construct a biounit from an assymetric unit. Can be serialized.

  :param mmcif_buinfo: Biounit definition
  :type mmcif_buinfo: :class:`ost.io.MMCifInfoBioUnit`

  .. method:: ToBytes()

    :returns: A byte string from which the object can be reconstructed.

  .. staticmethod:: FromBytes(byte_string)

    :param byte_string: Can be created with :func:`ToBytes`
    :returns: A :class:`BUInfo` object

.. function:: CreateBU(asu, bu_info)

  Constructs a biounit given an assymetric unit and transformation
  information. The following properties are copied from the assymetric
  unit and are expected to be set (this is the case if you use
  :func:`ost.io.LoadMMCIF` for the assymetric unit):

  * Chain level:

    * Chain type (see :attr:`ost.mol.ChainHandle.type`)

  * Residue level:

    * Chem type (see :attr:`ost.mol.ResidueHandle.chem_type`)
    * Chem class (:attr:`ost.mol.ResidueHandle.chem_class`)
    * One letter code (see :attr:`ost.mol.ResidueHandle.one_letter_code`)
    * Secondary structure (see :attr:`ost.mol.ResidueHandle.sec_structure`)
    * IsProtein and IsLigand properties (see :attr:`ost.mol.ResidueHandle.is_protein`/:attr:`ost.mol.ResidueHandle.is_ligand`)

  Each chain in the returned biounit can be referenced back to the
  assymetric unit as they follow a standardised naming scheme:
  <*idx*>.<*asu_cname*>, where *asu_cname* is the chain name in the assymetric
  unit and *idx* is the nth occurence of that chain in the biounit with 
  one based indexing. There is a quirk though. An index of 1, for example 1.A,
  is reserved for the original AU chain with identity transform (read: no
  transform) applied. If a certain AU chain only occurs with an actual
  transform applied, numbering starts at 2.
  
  .. warning::
    There is the (rare) possibility that a AU chain that has only identity
    transform applied is not named 1.<au_cname>.
    As of january 2024, there are 3 pdb entries (8qn6, 8x1h, 2c0x) where
    the same AU chain with identity transform occurs several times in the same
    biounit. This is likely an error in the respective mmCIF files as the
    resulting chains sit on top of each other. OST just names the FIRST
    occurence as 1.<au_cname>.
    

  :param asu: The assymetric unit
  :type asu: :class:`ost.mol.EntityHandle`
  :param bu_info: Info object
  :type bu_info: :class:`MMCifInfoBioUnit`/:class:`BUInfo`
  :returns: A :class:`ost.mol.EntityHandle` of the requested biounit 