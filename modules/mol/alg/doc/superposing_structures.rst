.. currentmodule:: ost.mol.alg

Superposing structures
================================================================================

.. autofunction:: Superpose

.. autofunction:: ParseAtomNames

.. autofunction:: MatchResidueByNum

.. autofunction:: MatchResidueByIdx

.. autofunction:: MatchResidueByLocalAln

.. autofunction:: MatchResidueByGlobalAln

.. class:: SuperpositionResult

  .. attribute:: rmsd

    RMSD of the superposed entities.

  .. attribute:: view1
                 view2

    Two :class:`~ost.mol.EntityView` used in superposition (not set if methods
    with :class:`~ost.geom.Vec3List` used).

  .. attribute:: transformation

    Transformation (:class:`~ost.geom.Mat4`) used to map :attr:`view1` onto
    :attr:`view2`.

  .. attribute:: fraction_superposed
                 rmsd_superposed_atoms
                 ncycles

    For iterative superposition (:func:`IterativeSuperposeSVD`): fraction and
    RMSD of atoms that were superposed with a distance below the given
    threshold and the number of iteration cycles performed.

.. method:: SuperposeSVD(view1, view2, apply_transform=True)
            SuperposeSVD(list1, list2)

  Superposition of two sets of atoms minimizing RMSD using a classic SVD based
  algorithm.

  Note that the atom positions in the view are taken blindly in the order in
  which the atoms appear.

  :param view1: View on the model entity
  :type view1:  :class:`~ost.mol.EntityView`
  :param view2: View on the reference entity
  :type view2:  :class:`~ost.mol.EntityView`
  :param list1: List of atom positions for model entity
  :type list1:  :class:`~ost.geom.Vec3List`
  :param list2: List of atom positions for reference entity
  :type list2:  :class:`~ost.geom.Vec3List`
  :param apply_transform: If True, the superposition transform is applied to
                          the (full!) entity handle linked to *view1*.
  :type apply_transform:  :class:`bool`

  :return: An instance of :class:`SuperpositionResult`.

.. method:: IterativeSuperposeSVD(view1, view2, max_iterations=5, \
                                  distance_threshold=3.0, apply_transform=True)
            IterativeSuperposeSVD(list1, list2, max_iterations=5, \
                                  distance_threshold=3.0)

  Iterative superposition of two sets of atoms. In each iteration cycle, we
  keep a fraction of atoms with distances below *distance_threshold* and get
  the superposition considering only those atoms.

  Note that the atom positions in the view are taken blindly in the order in
  which the atoms appear.

  :param view1: View on the model entity
  :type view1:  :class:`~ost.mol.EntityView`
  :param view2: View on the reference entity
  :type view2:  :class:`~ost.mol.EntityView`
  :param list1: List of atom positions for model entity
  :type list1:  :class:`~ost.geom.Vec3List`
  :param list2: List of atom positions for reference entity
  :type list2:  :class:`~ost.geom.Vec3List`
  :param max_iterations: Max. number of iterations to be performed
  :type max_iterations:  :class:`int`
  :param distance_threshold: Distance threshold defining superposed atoms
  :type distance_threshold:  :class:`float`
  :param apply_transform: If True, the superposition transform is applied to
                          the (full!) entity handle linked to *view1*.
  :type apply_transform:  :class:`bool`

  :return: An instance of :class:`SuperpositionResult`.

  :raises: Exception if atom counts do not match or if less than 3 atoms.

.. method:: CalculateRMSD(view1, view2, transformation=geom.Mat4())

  :return: RMSD of atom positions (taken blindly in the order in which the
           atoms appear) in the two given views.
  :rtype:  :class:`float`

  :param view1: View on the model entity
  :type view1:  :class:`~ost.mol.EntityView`
  :param view2: View on the reference entity
  :type view2:  :class:`~ost.mol.EntityView`
  :param transformation: Optional transformation to apply on each atom position
                         of *view1*.
  :type transformation:  :class:`~ost.geom.Mat4` 