:mod:`~ost.mol.alg` -- Algorithms for Structures
================================================================================

.. module:: ost.mol.alg
   :synopsis: Algorithms operating on molecular structures

Topics
--------------------------------------------------------------------------------

.. toctree::
  :maxdepth: 1

  biounits
  gdt
  lddt
  mapping_functions
  molck
  steric_clashes
  superposing_structures


Submodules
--------------------------------------------------------------------------------

.. toctree::
  :maxdepth: 1

  chain_mapping
  contact_score
  dockq
  helix_kinks
  ligand_scoring
  qsscore
  scoring
  scoring_base
  stereochemistry
  structure_analysis
  trajectory_analysis

Algorithms on Structures
--------------------------------------------------------------------------------

.. method:: Accessibility(ent, probe_radius=1.4, include_hydrogens=False,\
                          include_hetatm=False, include_water=False,\
                          oligo_mode=False, selection="", asa_abs="asaAbs",\
                          asa_rel="asaRel", asa_atom="asaAtom", \
                          algorithm = NACCESS)
            
  Calculates the accesssible surface area for ever atom in *ent*. The algorithm
  mimics the behaviour of the bindings available for the NACCESS and DSSP tools 
  and has been tested to reproduce the numbers accordingly.

  :param ent:           Entity on which to calculate the surface
  :type ent:            :class:`~ost.mol.EntityView` /
                        :class:`~ost.mol.EntityHandle`

  :param probe_radius:  Radius of probe to determine accessible surface area
  :type probe_radius:  :class:`float`

  :param include_hydrogens: Whether to include hydrogens in the solvent
                            accessibility calculations. By default, 
                            every atom with ele=H,D is simply neglected.
  :type include_hydrogens:  :class:`bool`

  :param include_hetatms: Whether to include atoms flagged as hetatms 
                          , i.e. ligands, in the solvent
                          accessibility calculations. They are neglected 
                          by default.
  :type include_hetatms:  :class:`bool`

  :param include_water: Whether to include water in the solvent
                        accessibility calculations. By default, 
                        every residue with name "HOH" is neglected.
  :type include_water:  :class:`bool`

  :param oligo_mode:    A typical used case of accessibility calculations is to 
                        determine the solvent accessibility of a full complex
                        and then the accessibility of each chain individually.
                        Lots of calculations can be cached because only the 
                        values of the atoms close to an interface change.
                        This is exactly what happens when you activate the 
                        oligo mode. It returns exactly the same value but adds,
                        additionally to the values estimated in full complex, 
                        the values from each individual chain as float 
                        properties on every residue and atom. Example for atom 
                        accessible surface if the according property name is set 
                        to "asaAtom": Accessibility in the full complex is 
                        stored as "asaAtom", the accessibility when only 
                        considering that particular chain is stored as 
                        "asaAtom_single_chain".
                        The other properties follow the same naming scheme.
  :type oligo_mode:     :class:`bool`

  :param selection:     Selection statement, that gets applied on *ent* before 
                        doing anything. Everything that is not selected is 
                        neglected. The default value of "" results in no
                        selection at all.
  :type selection:      :class:`str`

  :param asa_abs:       Float property name to assign the summed solvent
                        accessible surface from each atom to a residue.
  :type asa_abs:        :class:`str`

  :param asa_rel:       Float property name to assign the relative solvent 
                        accessibility to a residue. This is the absolute 
                        accessibility divided by the maximum solvent 
                        accessibility of that particular residue. 
                        This maximum solvent accessibility is dependent on
                        the chosen :class:`AccessibilityAlgorithm`.
                        Only residues of the 20 standarad amino acids 
                        can be handled.

                        * In case of the **NACCESS** algorithm you can expect
                          a value in range [0.0, 100.0] and a value of -99.9
                          for non standard residues.

                        * In case of the **DSSP** algorithm you can expect a
                          value in range [0.0, 1.0], no float property is assigned
                          in case of a non standard residue.

  :type asa_rel:       :class:`str`

  :param asa_atom:      Float property name to assign the solvent accessible 
                        area to each atom.
  :type asa_atom:       :class:`str`   

  :param algorithm:     Specifies the used algorithm for solvent accessibility
                        calculations

  :type algorithm:      :class:`AccessibilityAlgorithm`   

  

  :return: The summed solvent accessibilty of each atom in *ent*.

.. class:: AccessibilityAlgorithm

  The accessibility algorithm enum specifies the algorithm used by
  the respective tools. Available are:

    NACCESS, DSSP



.. method:: AssignSecStruct(ent)

  Assigns secondary structures to all residues based on hydrogen bond patterns
  as described by DSSP.

  :param ent:           Entity on which to assign secondary structures
  :type ent:            :class:`~ost.mol.EntityView` /
                        :class:`~ost.mol.EntityHandle`



.. class:: FindMemParam

  Result object for the membrane detection algorithm described below

  .. attribute:: axis

    initial search axis from which optimal membrane slab could be found

  .. attribute:: tilt_axis

    Axis around which we tilt the membrane starting from the initial axis

  .. attribute:: tilt

    Angle to tilt around tilt axis

  .. attribute:: angle

    After the tilt operation we perform a rotation around the initial axis
    with this angle to get the final membrane axis

  .. attribute:: membrane_axis

    The result of applying the tilt and rotation procedure described above.
    The membrane_axis is orthogonal to the membrane plane and has unit length.

  .. attribute:: pos

    Real number that describes the membrane center point. To get the actual
    position you can do: pos * membrane_axis

  .. attribute:: width

    Total width of the membrane in A

  .. attribute:: energy

    Pseudo energy of the implicit solvation model

  .. attribute:: membrane_asa

    Membrane accessible surface area

  .. attribute:: membrane_representation

    Dummy atoms that represent the membrane. This entity is only valid if
    the according flag has been set to True when calling FindMembrane.


.. method:: FindMembrane(ent, assign_membrane_representation=True, fast=False)

  Estimates the optimal membrane position of a protein by using an implicit 
  solvation model. The original algorithm and the used energy function are 
  described in: Lomize AL, Pogozheva ID, Lomize MA, Mosberg HI (2006) 
  Positioning of proteins in membranes: A computational approach.

  There are some modifications in this implementation and the procedure is
  as follows:

  * Initial axis are constructed that build the starting point for initial 
    parameter grid searches.

  * For every axis, the protein is rotated so that the axis builds the z-axis

    * In order to exclude internal hydrophilic pores, only the outermost atoms
      with respect the the z-axis enter an initial grid search
    * The width and position of the membrane is optimized for different 
      combinations of tilt and rotation angles (further described in
      :class:`FindMemParam`). The top 20 parametrizations 
      (only top parametrization if *fast* is True) are stored for further 
      processing.

  * The 20 best membrane parametrizations from the initial grid search 
    (only the best if *fast* is set to True) enter a final 
    minimization step using a Levenberg-Marquardt minimizer. 


  :param ent:           Entity of a transmembrane protein, you'll get weird 
                        results if this is not the case. The energy term
                        of the result is typically a good indicator whether
                        *ent* is an actual transmembrane protein. The 
                        following float properties will be set on the atoms:

                        * 'asaAtom' on all atoms that are selected with 
                          ent.Select('peptide=true and ele!=H') as a result
                          of envoking :meth:`Accessibility`.
                        * 'membrane_e' the contribution of the potentially 
                          membrane facing atoms to the energy function. 
                          
  :type ent:            :class:`ost.mol.EntityHandle` / :class:`ost.mol.EntityView`

  :param assign_membrane_representation: Whether to construct a membrane 
                                         representation using dummy atoms

  :type assign_membrane_representation: :class:`bool`

  :param fast:          If set to false, the 20 best results of the initial grid
                        search undergo a Levenberg-Marquardt minimization and
                        the parametrization with optimal minimized energy is 
                        returned. 
                        If set to yes, only the best result of the initial grid
                        search is selected and returned after 
                        Levenberg-Marquardt minimization.

  :returns:             The results object
  :rtype:               :class:`ost.mol.alg.FindMemParam`

.. _traj-analysis:

Trajectory Analysis
--------------------------------------------------------------------------------

This is a set of functions used for basic trajectory analysis such as extracting
positions, distances, angles and RMSDs. The organization is such that most
functions have their counterpart at the individual :class:`frame level
<ost.mol.CoordFrame>` so that they can also be called on one frame instead of
the whole trajectory.

All these functions have a "stride" argument that defaults to stride=1, which is
used to skip frames in the analysis.

Additional functions are available in the :mod:`ost.mol.alg.trajectory_analysis`
submodule.


.. function:: AnalyzeAngle(traj, atom1, atom2, atom3, stride=1)

  This function extracts the angle between three atoms from a trajectory 
  and returns it as a vector. The second atom is taken as being the central
  atom, so that the angle is between the vectors (atom1.pos-atom2.pos) and
  (atom3.pos-atom2.pos).
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The first :class:`~ost.mol.AtomHandle`.
  :param atom2: The second :class:`~ost.mol.AtomHandle`.
  :param atom3: The third :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.


.. function:: AnalyzeAromaticRingInteraction(traj, view_ring1, view_ring2, stride=1)

  This function is a crude analysis of aromatic ring interactions. For each frame in a trajectory, it calculates
  the minimal distance between the atoms in one view and the center of mass of the other
  and vice versa, and returns the minimum between these two minimal distances.
  Concretely, if the two views are the heavy atoms of two rings, then it returns the minimal
  center of mass - heavy atom distance betweent he two rings

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param view_ring1: First group of atoms
  :type view_ring1: :class:`~ost.mol.EntityView`.
  :param view_ring2: Second group of atoms
  :type view_ring2: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.


.. function:: AnalyzeAtomPos(traj, atom1, stride=1)

  This function extracts the position of an atom from a trajectory. It returns 
  a vector containing the position of the atom for each analyzed frame.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.
  
  
.. function:: AnalyzeCenterOfMassPos(traj, sele, stride=1)

  This function extracts the position of the center-of-mass of a selection 
  (:class:`~ost.mol.EntityView`) from a trajectory and returns it as a vector.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param sele: The selection from which the center of mass is computed
  :type sele: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.
  

.. function:: AnalyzeDihedralAngle(traj, atom1, atom2, atom3, atom4, stride=1)

  This function extracts the dihedral angle between four atoms from a trajectory 
  and returns it as a vector. The angle is between the planes containing the  
  first three and the last three atoms.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The first :class:`~ost.mol.AtomHandle`.
  :param atom2: The second :class:`~ost.mol.AtomHandle`.
  :param atom3: The third :class:`~ost.mol.AtomHandle`.
  :param atom4: The fourth :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.

.. function:: AnalyzeDistanceBetwAtoms(traj, atom1, atom2, stride=1)

  This function extracts the distance between two atoms from a trajectory 
  and returns it as a vector.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The first :class:`~ost.mol.AtomHandle`.
  :param atom2: The second :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.  


.. function:: AnalyzeDistanceBetwCenterOfMass(traj, sele1, sele2, stride=1)

  This function extracts the distance between the center-of-mass of two 
  selections (:class:`~ost.mol.EntityView`) from a trajectory and returns it as 
  a vector.

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param sele1: The selection from which the first center of mass is computed
  :type sele1: :class:`~ost.mol.EntityView`.
  :param sele2: The selection from which the second center of mass is computed
  :type sele2: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.


.. function:: AnalyzeMinDistance(traj, view1, view2, stride=1)

  This function extracts the minimal distance between two sets of atoms 
  (view1 and view2) for each frame in a trajectory and returns it as a vector.

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param view1: The first group of atoms
  :type view1: :class:`~ost.mol.EntityView`.
  :param view2: The second group of atoms
  :type view2: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.
     
.. function:: AnalyzeMinDistanceBetwCenterOfMassAndView(traj, view_cm, view_atoms, stride=1)

  This function extracts the minimal distance between a set of atoms 
  (view_atoms) and the center of mass of a second set of atoms (view_cm) 
  for each frame in a trajectory and returns it as a vector.

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param view_cm: The group of atoms from which the center of mass is taken
  :type view_cm: :class:`~ost.mol.EntityView`.
  :param view_atoms: The second group of atoms
  :type view_atoms: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.


.. function:: AnalyzeRMSD(traj, reference_view, sele_view, stride=1)

  This function extracts the rmsd between two :class:`~ost.mol.EntityView` and 
  returns it as a vector. The views don't have to be from the same entity. The 
  reference positions are taken directly from the reference_view, evaluated only 
  once. The positions from the sele_view are evaluated for each frame.
  If you want to compare to frame i of the trajectory t, first use 
  t.CopyFrame(i) for example:
 
  .. code-block:: python
     
    eh = io.LoadPDB(...)
    t = io.LoadCHARMMTraj(eh, ...)
    sele = eh.Select(...)
    t.CopyFrame(0)
    mol.alg.AnalyzeRMSD(t, sele, sele)

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param reference_view: The selection used as reference structure
  :type reference_view: :class:`~ost.mol.EntityView`.
  :param sele_view: The selection compared to the reference_view
  :type sele_view: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.


.. function:: SuperposeFrames(frames, sel, from=0, to=-1, ref=-1)

  This function superposes the frames of the given coord group and returns them
  as a new coord group.
  
  :param frames: The source coord group.
  :type frames: :class:`~ost.mol.CoordGroupHandle`
  :param sel: An entity view containing the selection of atoms to be used for     
    superposition. If set to an invalid view, all atoms in the coord group are 
    used.
  :type sel: :class:`ost.mol.EntityView`
  :param from: index of the first frame
  :param to: index of the last frame plus one. If set to -1, the value is set to 
     the number of frames in the coord group
  :param ref: The index of the reference frame to use for superposition. If set 
     to -1, the each frame is superposed to the previous frame.
     
  :returns: A newly created coord group containing the superposed frames.

.. function:: SuperposeFrames(frames, sel, ref_view, from=0, to=-1)
  :noindex:

  Same as SuperposeFrames above, but the superposition is done on a reference
  view and not on another frame of the trajectory.
  
  :param frames: The source coord group.
  :type frames: :class:`~ost.mol.CoordGroupHandle`
  :param sel: An entity view containing the selection of atoms of the frames to be used for     
    superposition.
  :type sel: :class:`ost.mol.EntityView`
  :param ref_view: The reference view on which the frames will be superposed. The number
    of atoms in this reference view should be equal to the number of atoms in sel.
  :type ref_view: :class:`ost.mol.EntityView`
  :param from: index of the first frame
  :param to: index of the last frame plus one. If set to -1, the value is set to 
     the number of frames in the coord group     
  
  :returns: A newly created coord group containing the superposed frames.
