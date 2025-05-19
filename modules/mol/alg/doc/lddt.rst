Local Distance Difference Test (LDDT)
================================================================================

.. note::

  This is a new implementation of LDDT, introduced in OpenStructure 2.4 with
  focus on supporting quaternary structure and compounds beyond the 20
  standard proteinogenic amino acids.
  The :doc:`previous LDDT code <lddt_deprecated>` that comes with
  `Mariani et al. <https://dx.doi.org/10.1093/bioinformatics/btt473>`_ is
  considered deprecated.

.. note::

  :class:`lddt.lDDTScorer` provides the raw Python API to compute LDDT but
  stereochemistry checks as described in
  `Mariani et al. <https://dx.doi.org/10.1093/bioinformatics/btt473>`_
  must be done seperately. You may want to check out the
  ``compare-structures`` action (:ref:`ost compare structures`) to
  compute LDDT with pre-processing and support for quaternary structures.


.. autoclass:: ost.mol.alg.lddt.lDDTScorer
  :members:

.. autoclass:: ost.mol.alg.lddt.SymmetrySettings
  :members:

.. autofunction:: ost.mol.alg.lddt.GetDefaultSymmetrySettings

.. autoclass:: ost.mol.alg.lddt.CustomCompound
  :members:

.. class:: lDDTSettings(radius=15, \
                        sequence_separation=0, \
                        cutoffs=(0.5, 1.0, 2.0, 4.0), \
                        label="locallddt")

  Object containing the settings used for LDDT calculations.

  :param radius: Sets :attr:`radius`.
  :param sequence_separation: Sets :attr:`sequence_separation`.
  :param cutoffs: Sets :attr:`cutoffs`.
  :param label: Sets :attr:`label`.

  .. attribute:: radius

    Distance inclusion radius.

    :type: :class:`float`

  .. attribute:: sequence_separation

    Sequence separation.

    :type: :class:`int`

  .. attribute:: cutoffs

    List of thresholds used to determine distance conservation.

    :type: :class:`list` of :class:`float`

  .. attribute:: label

    The base name for the ResidueHandle properties that store the local scores.

    :type: :class:`str`

  .. method:: PrintParameters()

    Print settings.

  .. method:: ToString()

    :return: String representation of the lDDTSettings object.
    :rtype:  :class:`str` 