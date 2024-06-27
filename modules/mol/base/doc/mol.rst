:mod:`~ost.mol` -- Molecular structures and surfaces
================================================================================

.. module:: ost.mol
   :synopsis: Contains classes and functions to deal with molecular structures
              and surfaces

The mol module implements data structures to work with molecular datasets. At its heart lie the :class:`EntityHandle` and :class:`EntityView` classes which represent molecular structures such as proteins, DNA, RNA and small molecules. There are also classes to deal with molecular surfaces.

.. toctree::
  
  entity
  editors
  query
  surface
  traj
  ../alg/molalg
  ../alg/chain_mapping
  ../alg/contact_score
  ../alg/dockq
  ../alg/helix_kinks
  ../alg/ligand_scoring
  ../alg/qsscore
  ../alg/scoring
  ../alg/stereochemistry
  ../alg/structure_analysis
  ../alg/trajectory_analysis