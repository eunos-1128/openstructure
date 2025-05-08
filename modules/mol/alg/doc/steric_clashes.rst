Steric Clashes
================================================================================

The following function detects steric clashes in atomic structures. Two atoms are clashing if their euclidian distance is smaller than a threshold value (minus a tolerance offset). 


.. function:: FilterClashes(entity, clashing_distances, always_remove_bb=False)

  This function filters out residues with non-bonded clashing atoms. If the
  clashing atom is a backbone atom, the complete residue is removed from the
  structure, if the atom is part of the sidechain, only the sidechain atoms are
  removed. This behavior is changed  by the *always_remove_bb* flag: when the
  flag is set to True the whole residue is removed even if a clash is just
  detected in the side-chain.

  The function returns a view containing all elements (residues, atoms) that
  have not been removed from the input structure, plus a
  :class:`~ost.mol.alg.ClashingInfo` object containing information about the
  detected clashes.
  
  Two atoms are defined as clashing if their distance is shorter than the
  reference distance minus a tolerance threshold. The information about the
  clashing distances and the tolerance thresholds for all possible pairs of
  atoms is passed to the function as a parameter.

  Hydrogen and deuterium atoms are ignored by this function.
  
  :param entity: The input entity
  :type entity: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param clashing_distances: information about the clashing distances
  :type clashing_distances: :class:`~ost.mol.alg.ClashingDistances`
  :param always_remove_bb: if set to True, the whole residue is removed even if
                           the clash happens in the side-chain
  :type always_remove_bb:  :class:`bool`

  :returns: A tuple of two elements: The filtered :class:`~ost.mol.EntityView`,
            and a :class:`~ost.mol.alg.ClashingInfo` object


.. function:: CheckStereoChemistry(entity, bond_stats, angle_stats, \
                                   bond_tolerance, angle_tolerance, \
                                   always_remove_bb=False)

  This function filters out residues with severe stereo-chemical violations. If
  the violation involves a backbone atom, the complete residue is removed from
  the structure, if it involves an atom that is part of the sidechain, only the
  sidechain is removed. This behavior is changed  by the *always_remove_bb*
  flag: when the flag is set to True the whole residue is removed even if a
  violation is just detected in the side-chain.

  The function returns a view containing all elements (residues, atoms) that
  have not been removed from the input structure, plus a
  :class:`~ost.mol.alg.StereoChemistryInfo` object containing information about
  the detected stereo-chemical violations.
    
  A violation is defined as a bond length that lies outside of the range:
  [mean_length-std_dev*bond_tolerance, mean_length+std_dev*bond_tolerance] or an
  angle width outside of the range [mean_width-std_dev*angle_tolerance,
  mean_width+std_dev*angle_tolerance ]. The information about the mean lengths
  and widths and the corresponding standard deviations is passed to the function
  using two parameters.

  Hydrogen and deuterium atoms are ignored by this function.

  :param entity: The input entity
  :type entity: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param bond_stats: statistics about bond lengths
  :type bond_stats: :class:`~ost.mol.alg.StereoChemicalParams`
  :param angle_stats: statistics about angle widths
  :type angle_stats: :class:`~ost.mol.alg.StereoChemicalParams`
  :param bond_tolerance: tolerance for bond lengths (in standard deviations)
  :type bond_tolerance:  :class:`float`
  :param angle_tolerance: tolerance for angle widths (in standard deviations)
  :type angle_tolerance:  :class:`float`
  :param always_remove_bb: if set to True, the whole residue is removed even if
                           the clash happens in the side-chain
  :type always_remove_bb:  :class:`bool`

  :returns: A tuple of two elements: The filtered :class:`~ost.mol.EntityView`, and a :class:`~ost.mol.alg.StereoChemistryInfo` object


.. class:: ClashingInfo

  This object is returned by the :func:`FilterClashes` function, and contains
  information about the clashes detected by the function.

  .. method:: GetClashCount()

    :return: number of clashes between non-bonded atoms detected in the
             input structure

  .. method:: GetAverageOffset()

    :return: a value in Angstroms representing the average offset by which
             clashing atoms lie closer than the minimum acceptable distance
             (which of course differs for each possible pair of elements)

  .. method:: GetClashList()

    :return: list of detected inter-atomic clashes
    :rtype:  :class:`list` of :class:`ClashEvent`


.. class:: ClashEvent

  This object contains all the information relative to a single clash detected
  by the :func:`FilterClashes` function

  .. method:: GetFirstAtom()
              GetSecondAtom()

    :return: atoms which clash
    :rtype:  :class:`~ost.mol.alg.UniqueAtomIdentifier`

  .. method:: GetModelDistance()

    :return: distance (in Angstroms) between the two clashing atoms as observed
             in the model

  .. method:: GetAdjustedReferenceDistance()

    :return: minimum acceptable distance (in Angstroms) between the two atoms
             involved in the clash, as defined in :class:`ClashingDistances`


.. class:: StereoChemistryInfo

  This object is returned by the :func:`CheckStereoChemistry` function, and
  contains information about bond lengths and planar angle widths in the
  structure that diverge from the parameters tabulated by Engh and Huber in the
  International Tables of Crystallography. Only elements that diverge from the
  tabulated value by a minimumnumber of standard deviations (defined when the
  CheckStereoChemistry function is called) are reported.

  .. method:: GetBadBondCount()

    :return: number of bonds where a serious violation was detected

  .. method:: GetBondCount()

    :return: total number of bonds in the structure checked by the
             CheckStereoChemistry function

  .. method:: GetAvgZscoreBonds()

    :return: average z-score of all the bond lengths in the structure, computed
             using Engh and Huber's mean and standard deviation values

  .. method:: GetBadAngleCount()

    :return: number of planar angles where a serious violation was detected

  .. method:: GetAngleCount()

    :return: total number of planar angles in the structure checked by the
             CheckStereoChemistry function

  .. method:: GetAvgZscoreAngles()

    :return: average z-score of all the planar angle widths, computed using Engh
             and Huber's mean and standard deviation values.

  .. method:: GetBondViolationList()

     :return: list of bond length violations detected in the structure
     :rtype:  :class:`list` of :class:`~ost.mol.alg.StereoChemicalBondViolation`

  .. method:: GetAngleViolationList()

     :return: list of angle width violations detected in the structure
     :rtype: :class:`list` of :class:`~ost.mol.alg.StereoChemicalAngleViolation`


.. class:: StereoChemicalBondViolation

  This object contains all the information relative to a single detected violation of stereo-chemical parameters in a bond length

  .. method:: GetFirstAtom()
              GetSecondAtom()

    :return: first / second atom of the bond
    :rtype:  :class:`~ost.mol.alg.UniqueAtomIdentifier`

  .. method:: GetBondLength()

    :return: length of the bond (in Angstroms) as observed in the model

  .. method:: GetAllowedRange()

    :return: allowed range of bond lengths (in Angstroms), according to Engh and
             Huber's tabulated parameters and the tolerance threshold used when
             the :func:`CheckStereoChemistry` function was called
    :rtype:  :class:`tuple` (minimum and maximum)


.. class:: StereoChemicalAngleViolation

  This object contains all the information relative to a single detected violation of stereo-chemical parameters in a planar angle width

  .. method:: GetFirstAtom()
              GetSecondAtom()
              GetThirdAtom()

    :return: first / second (vertex) / third atom that defines the planar angle
    :rtype:  :class:`UniqueAtomIdentifier`

  .. method:: GetAngleWidth()

    :return: width of the planar angle (in degrees) as observed in the model

  .. method:: GetAllowedRange()

    :return: allowed range of angle widths (in degrees), according to Engh and
             Huber's tabulated parameters and the tolerance threshold used when
             the :func:`CheckStereoChemistry` function was called
    :rtype:  :class:`tuple` (minimum and maximum)


.. class:: ClashingDistances

  Object containing information about clashing distances between non-bonded atoms

  .. method:: ClashingDistances()

    Creates an empty distance list

  .. method:: SetClashingDistance(ele1, ele2, clash_distance, tolerance)

    Adds or replaces an entry in the list

    :param ele1: string containing the first element's name 
    :param ele2: string containing the second element's name 
    :param clash_distance: minimum clashing distance (in Angstroms)
    :param tolerance: tolerance threshold (in Angstroms)

  .. method:: GetClashingDistance(ele1, ele2)

    :return: reference distance and a tolerance threshold (both in Angstroms)
             for two elements
    :rtype:  :class:`tuple` (minimum clashing distance, tolerance threshold)
    :param ele1: string containing the first element's name 
    :param ele2: string containing the second element's name 

  .. method:: GetAdjustedClashingDistance(ele1, ele2)

    :return: reference distance (in Angstroms) for two elements, already
             adjusted by the tolerance threshold
    :param ele1: string containing the first element's name
    :param ele2: string containing the second element's name

  .. method:: GetMaxAdjustedDistance()
 
    :return: longest clashing distance (in Angstroms) in the list, after
             adjustment with tolerance threshold

  .. method:: IsEmpty()

    :return: True if the list is empty (i.e. in an invalid, useless state)
 
  .. method:: PrintAllDistances()

    Prints all distances in the list to standard output


.. class:: StereoChemicalParams

  Object containing stereo-chemical information about bonds and angles. For each
  item (bond or angle in a specific residue), stores the mean and standard
  deviation

  .. method:: StereoChemicalParams()

    Creates an empty parameter list

  .. method:: SetParam(item, residue, mean, standard_dev)

    Adds or replaces an entry in the list

    :param item: string defining a bond (format: X-Y) or an angle (format:
                 X-Y-Z), where X,Y an Z are atom names
    :param residue: string containing the residue type for this entry
    :param mean: mean bond length (in Angstroms) or angle width (in degrees)
    :param standard_dev: standard deviation of the bond length (in Angstroms)
                         or of the angle width (in degrees)

  .. method GetParam(item, residue)

    :return: entry from the list as set in :meth:`SetParam`
    :rtype:  :class:`tuple` (mean, standard deviation)
    :param item: string as used in :meth:`SetParam`
    :param residue: string as used in :meth:`SetParam`

  .. method ContainsParam(item, residue)

    :return: True if a specific entry is present in the list, False if not
    :param item: string as used in :meth:`SetParam`
    :param residue: string as used in :meth:`SetParam`

  .. method:: IsEmpty()

    :return: True if the list is empty (i.e. in an invalid, useless state)
 
  .. method:: PrintAllParameters()

    Prints all entries in the list to standard output  


.. function:: FillClashingDistances(file_content)
              FillBondStereoChemicalParams(file_content)
              FillAngleStereoChemicalParams(file_content)

  These three functions fill a list of reference clashing distances, a list of
  stereo-chemical parameters for bonds and a list of stereo-chemical parameters
  for angles, respectively, starting from the content of a parameter file.

  :param file_content: list of lines from the parameter file
  :type file_content:  :class:`list` of :class:`str`

  :rtype: :class:`~ost.mol.alg.ClashingDistances` or
          :class:`~ost.mol.alg.StereoChemicalParams`


.. function:: FillClashingDistancesFromFile(filename)
              FillBondStereoChemicalParamsFromFile(filename)
              FillAngleStereoChemicalParamsFromFile(filename)

  These three functions fill a list of reference clashing distances, a list of
  stereo-chemical parameters for bonds and a list of stereo-chemical parameters
  for angles, respectively, starting from a file path.

  :param filename: path to parameter file
  :type filename:  :class:`str`

  :rtype: :class:`~ost.mol.alg.ClashingDistances` or
          :class:`~ost.mol.alg.StereoChemicalParams`


.. function:: DefaultClashingDistances()
              DefaultBondStereoChemicalParams()
              DefaultAngleStereoChemicalParams()

  These three functions fill a list of reference clashing distances, a list of
  stereo-chemical parameters for bonds and a list of stereo-chemical parameters
  for angles, respectively, using the default parameter files distributed with
  OpenStructure.

  :rtype: :class:`~ost.mol.alg.ClashingDistances` or
          :class:`~ost.mol.alg.StereoChemicalParams` 