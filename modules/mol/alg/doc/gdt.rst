GDT - Global Distance Test
================================================================================

Implements the GDT score, i.e. identifies the largest number of positions
that can be superposed within a given distance threshold. The final
GDT score is then the returned number divided by the total number of
reference positioons. The algorithm is similar to what is described for
the LGA tool but simpler. Therefore, the fractions reported by OpenStructure
tend to be systematically lower. For benchmarking we computed the full GDT_TS,
i.e. average GDT for distance thresholds [1, 2, 4, 8], on all CASP15 TS
models. 96.5% of differences to the LGA results from the predictioncenter are
within 2 GDT points and 99.2% are within 3 GDT points. The max difference
is 7.39 GDT points.

The algorithm expects two position lists of same length and applies a sliding
window with specified length to define a subset of position pairs as starting
point for iterative superposition. Each iterative superposition applies the
following steps:

- Compute minimal RMSD superposition on subset of position pairs
- Apply superposition on all model positions
- Compute pairwise distances of all model positions and reference positions
- Define new subset of position pairs: pairs within distance threshold
- Stop if subset doesn't change anymore

The subset in any of the iterations which is largest is stored.

This is done for each sliding window position and the largest subset ever
observed is reported. To avoid long runtimes for large problem sizes, the
sliding window is not applied on each possible position but is capped.
If the number of positions is larger than this threshold, the sliding
window is only applied on N equidistant locations.

.. function:: GDT(mdl_pos, ref_pos, window_size, max_windows, distance_thresh)
  
    Returns number of positions that can be superposed within
    *distance_thresh* and the respective transformation matrix.

    :param mdl_pos: Positions representing the model, typically alpha-carbon
                    positions
    :param ref_pos: Positions representing the reference, typically
                    alpha-carbon positions
    :param window_size: Size of the sliding window that is used to serve as
                        starting point for iterative superposition.
                        The described benchmark was done with a value of 7.
    :param max_windows: Cap for number of starting points. The described
                        benchmark was done with a value of 1000.
    :param distance_thresh: Distance threshold for GDT algorithm
    :type mdl_pos: :class:`ost.geom.Vec3List`
    :type ref_pos: :class:`ost.geom.Vec3List`
    :type window_size: :class:`int`
    :type max_windows: :class:`int`
    :type distance_thresh: :class:`float`
    :returns: :class:`tuple` with first element being the number of
              superposable positions (:class:`int`) and the second element the
              transformation matrix (:class:`ost.geom.Mat4`) 