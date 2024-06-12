Parasail - Pairwise sequence alignments
================================================================================

Basic access to pairwise sequence alignment functionality in
`parasail <https://github.com/jeffdaily/parasail/>`_ from
within OpenStructure.

Citation:

  Jeff Daily. Parasail: SIMD C library for global, semi-global,
  and local pairwise sequence alignments. (2016) BMC Bioinformatics

Parasail support must be enabled at compile time - see installation
instructions. Parasail allows to choose from various strategies but
for the sake of simplicity, this Python binding always calls
``parasail_<mode>_trace_scan_sat`` which seems reasonably fast
across the global, semi-global and local modes.
See parasail documentation for more information. Example usage:

.. code-block:: python

  from ost import bindings

  print("parasail available?", bindings.ParasailAvailable())

  s1 = seq.CreateSequence("A", "HEAGAWGHEE")
  s2 = seq.CreateSequence("B", "PAWHEA")

  aln = bindings.ParaGlobalAlign(s1, s2, gap_open_penalty=11,
                                 gap_extension_penalty=1,
                                 matrix="blosum62")
                               
  print(aln.ToString())


.. method:: ParasailAvailable()

  returns if OpenStructure has been compiled with parasail support

.. method:: ParaGlobalAlign(s1, s2, gap_open_penalty=11, \
                            gap_extension_penalty=1, matrix="blosum62")

  Calls parasail_nw_trace_scan_sat to perform a Needleman-Wunsch global
  alignment.

  :param s1: First sequence
  :type s1: :class:`ost.seq.SequenceHandle`
  :param s2: Second sequence
  :type s2: :class:`ost.seq.SequenceHandle`
  :param gap_open_penalty: Gap opening penalty for dynamic programming table
  :type gap_open_penalty: :class:`int`
  :param gap_extension_penalty: Gap extension penalty for dynamic programming
                                table
  :type gap_extension_penalty: :class:`int`
  :param matrix: Substution matrix - Anything that can be resolved by parasail.
                 E.g. pretty much the full BLOSUM/PAM families of matrices
                 ("blosum62", "pam100",...) for protein sequences, or
                 "nuc44"/"dnafull" for nucleotides.
  :returns: :class:`ost.seq.AlignmentHandle`


.. method:: ParaSemiGlobalAlign(s1, s2, gap_open_penalty=11, \
                                gap_extension_penalty=1, matrix="blosum62")

  Same as :meth:`ParaGlobalAlign` but calling parasail_sg_trace_scan_sat to
  perform a semi-global alignment.


.. method:: ParaLocalAlign(s1, s2, gap_open_penalty=11, \
                            gap_extension_penalty=1, matrix="blosum62")

  Same as :meth:`ParaGlobalAlign` but calling parasail_sw_trace_scan_sat to
  perform a Smith-Waterman local alignment. The sequences in the returned
  alignment might be subsequences of *s1*/*s2* but have offsets set.
