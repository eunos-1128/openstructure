:mod:`~ost.seq.alg.aaindex` -- AAIndex annotations
--------------------------------------------------------------------------------

.. module:: ost.seq.alg.aaindex
   :synopsis: AAIndex annotations

.. autoclass:: ost.seq.alg.aaindex.AAIndex
  :members:
  :special-members: __getitem__

The annotations/scores can either refer to single amino acids or represent
pairwise values. The two types are:

.. autoclass:: ost.seq.alg.aaindex.AnnoType
  :members:
  :undoc-members:

The actual data of an entry in the aaindex database is stored in a
:class:`aaindex.AAIndexData` object:

.. autoclass:: ost.seq.alg.aaindex.AAIndexData
  :members:
