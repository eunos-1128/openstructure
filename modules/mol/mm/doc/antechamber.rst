:mod:`~ost.mol.mm.antechamber` -- Generating forcefields with Antechamber
--------------------------------------------------------------------------------

The antechamber submodule of mol.mm defines functions to use Antechamber (from
AmberTools) to automatically generate force field parameters and load the
results into :class:`~ost.mol.mm.Forcefield` objects.

**Example usage**:

.. code-block:: python

  from ost.mol import mm

  # create parameters for RVP using PDB's component dictionary
  mm.antechamber.RunAntechamber('RVP', 'components.cif', base_out_dir='ligands')

  # create force field
  ff = mm.Forcefield()
  ff = mm.antechamber.AddFromPath(ff, 'ligands/RVP')
  # equivalent: ff = mm.antechamber.AddFromFiles(ff, 'ligands/RVP/frcmod',
  #                                              'ligands/RVP/out.mpdb')
  # since Antechamber cannot deal with ions, you can do it manually
  ff = mm.antechamber.AddIon(ff, 'CL', 'CL', 35.45, -1.0, 0.4401, 0.4184)
  # save it
  ff.Save('ligands/ff.dat')

**Functions**:

.. automodule:: ost.mol.mm.antechamber
  :members: