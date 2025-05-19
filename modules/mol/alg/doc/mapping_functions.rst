Mapping Functions
================================================================================

These functions convert one residue into another. They are primarily used for standard and modified amino acids.

.. function:: CopyResidue(src_res, dst_res, editor)

  This function copies atoms from a source residue to a destination residue. It calls either
  :func:`CopyConserved` or :func:`CopyNonConserved` depending on whether the residues are
  identical or not.

  :param src_res: The source residue.
  :type src_res: :class:`~ost.mol.ResidueHandle`
  :param dst_res: The destination residue.
  :type dst_res: :class:`~ost.mol.ResidueHandle`
  :param editor: The editor to use for the copy operation.
  :type editor: :class:`~ost.mol.EntityEditor`

  :returns: The number of atoms copied.


.. function:: CopyConserved(src_res, dst_res, editor)

  This function copies atoms under the assumption that the source and destination residues
  are identical or that the source is a modified version of the destination. It handles
  modifications and converts selenium to sulfur for specific residues.

  :param src_res: The source residue.
  :type src_res: :class:`~ost.mol.ResidueHandle`
  :param dst_res: The destination residue.
  :type dst_res: :class:`~ost.mol.ResidueHandle`
  :param editor: The editor to use for the copy operation.
  :type editor: :class:`~ost.mol.EntityEditor`

  :returns: The number of atoms copied.


.. function:: CopyNonConserved(src_res, dst_res, editor)

  This function copies heavy backbone atoms and CBeta from the source to the destination.

  :param src_res: The source residue.
  :type src_res: :class:`~ost.mol.ResidueHandle`
  :param dst_res: The destination residue.
  :type dst_res: :class:`~ost.mol.ResidueHandle`
  :param editor: The editor to use for the copy operation.
  :type editor: :class:`~ost.mol.EntityEditor`

  :returns: The number of atoms copied. 