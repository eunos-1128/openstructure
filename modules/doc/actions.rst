..  Note on large code blocks: keep max. width to 100 or it will look bad
                               on webpage!
..  TODO: look at argparse directive to autogenerate --help output!

.. ost-actions:

OST Actions
================================================================================

A pure command line interface of OST is provided by actions.
You can execute ``ost -h`` for a list of possible actions and for every action,
you can type ``ost <ACTION> -h`` to get a description on its usage.

Here we list the most prominent actions with simple examples.

.. _ost compare structures:

Comparing two structures
--------------------------------------------------------------------------------

You can compare two structures from the command line with the
``ost compare-structures`` action. This can be considered a command line
interface to the :class:`~ost.mol.alg.scoring.Scorer`.

.. note::

  This is a new implementation of the ``compare-structures`` action, introduced
  in OpenStructure 2.4. The old version is still available as
  :doc:`compare-structures-legacy <deprecated_actions>`.

Details on the usage (output of ``ost compare-structures --help``):

.. code-block:: console

  usage: ost compare-structures [-h] -m MODEL -r REFERENCE [-o OUTPUT]
                                [-mf {pdb,cif,mmcif}] [-rf {pdb,cif,mmcif}]
                                [-mb MODEL_BIOUNIT] [-rb REFERENCE_BIOUNIT]
                                [-rna] [-ec] [-d] [-ds DUMP_SUFFIX] [-ft]
                                [-c CHAIN_MAPPING [CHAIN_MAPPING ...]] [--lddt]
                                [--local-lddt] [--bb-lddt] [--bb-local-lddt]
                                [--ilddt] [--cad-score] [--local-cad-score]
                                [--cad-exec CAD_EXEC]
                                [--usalign-exec USALIGN_EXEC]
                                [--override-usalign-mapping] [--qs-score]
                                [--dockq] [--dockq-capri-peptide] [--ics]
                                [--ips] [--rigid-scores] [--patch-scores]
                                [--tm-score] [--lddt-no-stereochecks]
                                [--n-max-naive N_MAX_NAIVE]
                                [--dump-aligned-residues] [--dump-pepnuc-alns]
                                [--dump-pepnuc-aligned-residues]
                                [--min-pep-length MIN_PEP_LENGTH]
                                [--min-nuc-length MIN_NUC_LENGTH] [-v VERBOSITY]
                                [--lddt-add-mdl-contacts]
  
  Evaluate model against reference 
  
  Example: ost compare-structures -m model.pdb -r reference.cif
  
  Loads the structures and performs basic cleanup:
  
   * Assign elements according to the PDB Chemical Component Dictionary
   * Map nonstandard residues to their parent residues as defined by the PDB
     Chemical Component Dictionary, e.g. phospho-serine => serine
   * Remove hydrogens
   * Remove OXT atoms
   * Remove unknown atoms, i.e. atoms that are not expected according to the PDB
     Chemical Component Dictionary
   * Select for peptide/nucleotide residues
  
  The cleaned structures are optionally dumped using -d/--dump-structures
  
  Output is written in JSON format (default: out.json). In case of no additional
  options, this is a dictionary with 8 keys describing model/reference comparison:
  
   * "reference_chains": Chain names of reference
   * "model_chains": Chain names of model
   * "chem_groups": Groups of polypeptides/polynucleotides from reference that
     are considered chemically equivalent. You can derive stoichiometry from this.
     Contains only chains that are considered in chain mapping, i.e. pass a
     size threshold (defaults: 6 for peptides, 4 for nucleotides).
   * "chem_mapping": List of same length as "chem_groups". Assigns model chains to
     the respective chem group. Again, only contains chains that are considered
     in chain mapping.
   * "chain_mapping": A dictionary with reference chain names as keys and the
     mapped model chain names as values. Missing chains are either not mapped
     (but present in "chem_groups", "chem_mapping") or were not considered in
     chain mapping (short peptides etc.)
   * "aln": Pairwise sequence alignment for each pair of mapped chains in fasta
     format.
   * "inconsistent_residues": List of strings that represent name mismatches of
     aligned residues in form
     <trg_cname>.<trg_rnum>.<trg_ins_code>-<mdl_cname>.<mdl_rnum>.<mdl_ins_code>.
     Inconsistencies may lead to corrupt results but do not abort the program.
     Program abortion in these cases can be enforced with
     -ec/--enforce-consistency.
   * "status": SUCCESS if everything ran through. In case of failure, the only
     content of the JSON output will be "status" set to FAILURE and an
     additional key: "traceback".
  
  The following additional keys store relevant input parameters to reproduce
  results:
  
   * "model"
   * "reference"
   * "fault_tolerant"
   * "model_biounit"
   * "reference_biounit"
   * "residue_number_alignment"
   * "enforce_consistency"
   * "cad_exec"
   * "usalign_exec"
   * "lddt_no_stereochecks"
   * "min_pep_length"
   * "min_nuc_length"
   * "lddt_add_mdl_contacts"
   * "dockq_capri_peptide"
  
  The pairwise sequence alignments are computed with Needleman-Wunsch using
  BLOSUM62 (NUC44 for nucleotides). Many benchmarking scenarios preprocess the
  structures to ensure matching residue numbers (CASP/CAMEO). In these cases,
  enabling -rna/--residue-number-alignment is recommended.
  
  Each score is opt-in and can be enabled with optional arguments.
  
  Example to compute global and per-residue lDDT values as well as QS-score:
  
  ost compare-structures -m model.pdb -r reference.cif --lddt --local-lddt --qs-score
  
  Example to inject custom chain mapping
  
  ost compare-structures -m model.pdb -r reference.cif -c A:B B:A
  
  options:
    -h, --help            show this help message and exit
    -m MODEL, --model MODEL
                          Path to model file.
    -r REFERENCE, --reference REFERENCE
                          Path to reference file.
    -o OUTPUT, --output OUTPUT
                          Output file name. The output will be saved as a JSON
                          file. default: out.json
    -mf {pdb,cif,mmcif}, --model-format {pdb,cif,mmcif}
                          Format of model file. pdb reads pdb but also pdb.gz,
                          same applies to cif/mmcif. Inferred from filepath if
                          not given.
    -rf {pdb,cif,mmcif}, --reference-format {pdb,cif,mmcif}
                          Format of reference file. pdb reads pdb but also
                          pdb.gz, same applies to cif/mmcif. Inferred from
                          filepath if not given.
    -mb MODEL_BIOUNIT, --model-biounit MODEL_BIOUNIT
                          Only has an effect if model is in mmcif format. By
                          default, the asymmetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the ID (as a string) of the one which
                          should be used.
    -rb REFERENCE_BIOUNIT, --reference-biounit REFERENCE_BIOUNIT
                          Only has an effect if reference is in mmcif format. By
                          default, the asymmetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the ID (as a string) of the one which
                          should be used.
    -rna, --residue-number-alignment
                          Make alignment based on residue number instead of
                          using a global BLOSUM62-based alignment (NUC44 for
                          nucleotides).
    -ec, --enforce-consistency
                          Enforce consistency. By default residue name
                          discrepancies between a model and reference are
                          reported but the program proceeds. If this flag is ON,
                          the program fails for these cases.
    -d, --dump-structures
                          Dump cleaned structures used to calculate all the
                          scores as PDB or mmCIF files using specified suffix.
                          Files will be dumped to the same location and in the
                          same format as original files.
    -ds DUMP_SUFFIX, --dump-suffix DUMP_SUFFIX
                          Use this suffix to dump structures. Defaults to
                          _compare_structures
    -ft, --fault-tolerant
                          Fault tolerant parsing.
    -c CHAIN_MAPPING [CHAIN_MAPPING ...], --chain-mapping CHAIN_MAPPING [CHAIN_MAPPING ...]
                          Custom mapping of chains between the reference and the
                          model. Each separate mapping consist of key:value
                          pairs where key is the chain name in reference and
                          value is the chain name in model.
    --lddt                Compute global lDDT score with default
                          parameterization and store as key "lddt".
                          Stereochemical irregularities affecting lDDT are
                          reported as keys "model_clashes", "model_bad_bonds",
                          "model_bad_angles" and the respective reference
                          counterparts.
    --local-lddt          Compute per-residue lDDT scores with default
                          parameterization and store as key "local_lddt". Score
                          for each residue is accessible by key
                          <chain_name>.<resnum>.<resnum_inscode>. Residue with
                          number 42 in chain X can be extracted with:
                          data["local_lddt"]["X.42."]. If there is an insertion
                          code, lets say A, the residue key becomes "X.42.A".
                          Stereochemical irregularities affecting lDDT are
                          reported as keys "model_clashes", "model_bad_bonds",
                          "model_bad_angles" and the respective reference
                          counterparts. Atoms specified in there follow the
                          following format:
                          <chain_name>.<resnum>.<resnum_inscode>.<atom_name>
    --bb-lddt             Compute global lDDT score with default
                          parameterization and store as key "bb_lddt". lDDT in
                          this case is only computed on backbone atoms: CA for
                          peptides and C3' for nucleotides
    --bb-local-lddt       Compute per-residue lDDT scores with default
                          parameterization and store as key "bb_local_lddt".
                          lDDT in this case is only computed on backbone atoms:
                          CA for peptides and C3' for nucleotides. Per-residue
                          scores are accessible as described for local_lddt.
    --ilddt               Compute global lDDT score which is solely based on
                          inter-chain contacts and store as key "ilddt". Same
                          stereochemical irregularities as for lddt apply.
    --cad-score           Compute global CAD's atom-atom (AA) score and store as
                          key "cad_score". --residue-number-alignment must be
                          enabled to compute this score. Requires
                          voronota_cadscore executable in PATH. Alternatively
                          you can set cad-exec.
    --local-cad-score     Compute local CAD's atom-atom (AA) scores and store as
                          key "local_cad_score". Per-residue scores are
                          accessible as described for local_lddt. --residue-
                          number-alignments must be enabled to compute this
                          score. Requires voronota_cadscore executable in PATH.
                          Alternatively you can set cad-exec.
    --cad-exec CAD_EXEC   Path to voronota-cadscore executable (installed from
                          https://github.com/kliment-olechnovic/voronota).
                          Searches PATH if not set.
    --usalign-exec USALIGN_EXEC
                          Path to USalign executable to compute TM-score. If not
                          given, an OpenStructure internal copy of USalign code
                          is used.
    --override-usalign-mapping
                          Override USalign mapping and inject our own rigid
                          mapping. Only works if external usalign executable is
                          provided that is reasonably new and contains that
                          feature.
    --qs-score            Compute QS-score, stored as key "qs_global", and the
                          QS-best variant, stored as key "qs_best". Interfaces
                          in the reference with non-zero contribution to QS-
                          score are available as key "qs_reference_interfaces",
                          the ones from the model as key "qs_model_interfaces".
                          "qs_interfaces" is a subset of
                          "qs_reference_interfaces" that contains interfaces
                          that can be mapped to the model. They are stored as
                          lists in format [ref_ch1, ref_ch2, mdl_ch1, mdl_ch2].
                          The respective per-interface scores for
                          "qs_interfaces" are available as keys
                          "per_interface_qs_global" and "per_interface_qs_best"
    --dockq               Compute DockQ scores and its components. Relevant
                          interfaces with at least one contact (any atom within
                          5A) of the reference structure are available as key
                          "dockq_reference_interfaces". Only interfaces between
                          peptide chains are considered here! Key
                          "dockq_interfaces" is a subset of
                          "dockq_reference_interfaces" that contains interfaces
                          that can be mapped to the model. They are stored as
                          lists in format [ref_ch1, ref_ch2, mdl_ch1, mdl_ch2].
                          The respective DockQ scores for "dockq_interfaces" are
                          available as key "dockq". It's components are
                          available as keys: "fnat" (fraction of reference
                          contacts which are also there in model) "irmsd"
                          (interface RMSD), "lrmsd" (ligand RMSD). The DockQ
                          score is strictly designed to score each interface
                          individually. We also provide two averaged versions to
                          get one full model score: "dockq_ave", "dockq_wave".
                          The first is simply the average of "dockq_scores", the
                          latter is a weighted average with weights derived from
                          number of contacts in the reference interfaces. These
                          two scores only consider interfaces that are present
                          in both, the model and the reference. "dockq_ave_full"
                          and "dockq_wave_full" add zeros in the average
                          computation for each interface that is only present in
                          the reference but not in the model.
    --dockq-capri-peptide
                          Flag that changes two things in the way DockQ and its
                          underlying scores are computed which is proposed by
                          the CAPRI community when scoring peptides (PMID:
                          31886916). ONE: Two residues are considered in contact
                          if any of their atoms is within 5A. This is relevant
                          for fnat and fnonat scores. CAPRI suggests to lower
                          this threshold to 4A for protein-peptide interactions.
                          TWO: irmsd is computed on interface residues. A
                          residue is defined as interface residue if any of its
                          atoms is within 10A of another chain. CAPRI suggests
                          to lower the default of 10A to 8A in combination with
                          only considering CB atoms for protein-peptide
                          interactions. Note that the resulting DockQ is not
                          evaluated for these slightly updated fnat and irmsd
                          (lrmsd stays the same).This flag has no influence on
                          patch_dockq scores.
    --ics                 Computes interface contact similarity (ICS) related
                          scores. A contact between two residues of different
                          chains is defined as having at least one heavy atom
                          within 5A. Contacts in reference structure are
                          available as key "reference_contacts". Each contact
                          specifies the interacting residues in format
                          "<cname>.<rnum>.<ins_code>". Model contacts are
                          available as key "model_contacts". The precision which
                          is available as key "ics_precision" reports the
                          fraction of model contacts that are also present in
                          the reference. The recall which is available as key
                          "ics_recall" reports the fraction of reference
                          contacts that are correctly reproduced in the model.
                          The ICS score (Interface Contact Similarity) available
                          as key "ics" combines precision and recall using the
                          F1-measure. All these measures are also available on a
                          per-interface basis for each interface in the
                          reference structure that are defined as chain pairs
                          with at least one contact (available as key
                          "contact_reference_interfaces"). The respective
                          metrics are available as keys
                          "per_interface_ics_precision",
                          "per_interface_ics_recall" and "per_interface_ics".
    --ips                 Computes interface patch similarity (IPS) related
                          scores. They focus on interface residues. They are
                          defined as having at least one contact to a residue
                          from any other chain. In short: if they show up in the
                          contact lists used to compute ICS. If ips is enabled,
                          these contacts get reported too and are available as
                          keys "reference_contacts" and "model_contacts".The
                          precision which is available as key "ips_precision"
                          reports the fraction of model interface residues, that
                          are also interface residues in the reference. The
                          recall which is available as key "ips_recall" reports
                          the fraction of reference interface residues that are
                          also interface residues in the model. The IPS score
                          (Interface Patch Similarity) available as key "ips" is
                          the Jaccard coefficient between interface residues in
                          reference and model. All these measures are also
                          available on a per-interface basis for each interface
                          in the reference structure that are defined as chain
                          pairs with at least one contact (available as key
                          "contact_reference_interfaces"). The respective
                          metrics are available as keys
                          "per_interface_ips_precision",
                          "per_interface_ips_recall" and "per_interface_ips".
    --rigid-scores        Computes rigid superposition based scores. They're
                          based on a Kabsch superposition of all mapped CA
                          positions (C3' for nucleotides). Makes the following
                          keys available: "oligo_gdtts": GDT with distance
                          thresholds [1.0, 2.0, 4.0, 8.0] given these positions
                          and transformation, "oligo_gdtha": same with
                          thresholds [0.5, 1.0, 2.0, 4.0], "rmsd": RMSD given
                          these positions and transformation, "transform": the
                          used 4x4 transformation matrix that superposes model
                          onto reference.
    --patch-scores        Local interface quality score used in CASP15. Scores
                          each model residue that is considered in the interface
                          (CB pos within 8A of any CB pos from another chain (CA
                          for GLY)). The local neighborhood gets represented by
                          "interface patches" which are scored with QS-score and
                          DockQ. Scores where not the full patches are
                          represented by the reference are set to None. Model
                          interface residues are available as key
                          "model_interface_residues", reference interface
                          residues as key "reference_interface_residues".
                          Residues are represented as string in form
                          <chain_name>.<resnum>.<resnum_inscode>. The respective
                          scores are available as keys "patch_qs" and
                          "patch_dockq"
    --tm-score            Computes TM-score with the USalign tool. Also computes
                          a chain mapping in case of complexes that is stored in
                          the same format as the default mapping. TM-score and
                          the mapping are available as keys "tm_score" and
                          "usalign_mapping"
    --lddt-no-stereochecks
                          Disable stereochecks for lDDT computation
    --n-max-naive N_MAX_NAIVE
                          Parameter for chain mapping. If the number of possible
                          mappings is <= *n_max_naive*, the full mapping
                          solution space is enumerated to find the the mapping
                          with optimal QS-score. A heuristic is used otherwise.
                          The default of 40320 corresponds to an octamer (8! =
                          40320). A structure with stoichiometry A6B2 would be
                          6!*2! = 1440 etc.
    --dump-aligned-residues
                          Dump additional info on aligned model and reference
                          residues.
    --dump-pepnuc-alns    Dump alignments of mapped chains but with sequences
                          that did not undergo Molck preprocessing in the
                          scorer. Sequences are extracted from model/target
                          after undergoing selection for peptide and nucleotide
                          residues.
    --dump-pepnuc-aligned-residues
                          Dump additional info on model and reference residues
                          that occur in pepnuc alignments.
    --min-pep-length MIN_PEP_LENGTH
                          Default: 6 - Relevant parameter if short peptides are
                          involved in scoring. Minimum peptide length for a
                          chain in the target structure to be considered in
                          chain mapping. The chain mapping algorithm first
                          performs an all vs. all pairwise sequence alignment to
                          identify "equal" chains within the target structure.
                          We go for simple sequence identity there. Short
                          sequences can be problematic as they may produce high
                          sequence identity alignments by pure chance.
    --min-nuc-length MIN_NUC_LENGTH
                          Default: 4 - Relevant parameter if short nucleotides
                          are involved in scoring.Minimum nucleotide length for
                          a chain in the target structure to be considered in
                          chain mapping. The chain mapping algorithm first
                          performs an all vs. all pairwise sequence alignment to
                          identify "equal" chains within the target structure.
                          We go for simple sequence identity there. Short
                          sequences can be problematic as they may produce high
                          sequence identity alignments by pure chance.
    -v VERBOSITY, --verbosity VERBOSITY
                          Set verbosity level. Defaults to 3 (Script).
    --lddt-add-mdl-contacts
                          Only using contacts in lDDT thatare within a certain
                          distance threshold in the reference does not penalize
                          for added model contacts. If set to True, this flag
                          will also consider reference contacts that are within
                          the specified distance threshold in the model but not
                          necessarily in the reference. No contact will be added
                          if the respective atom pair is not resolved in the
                          reference.

.. _ost compare ligand structures:

Comparing two structures with ligands
--------------------------------------------------------------------------------

You can compare two structures with non-polymer/small molecule ligands and
compute lDDT-PLI and ligand RMSD scores from the command line with the
``ost compare-ligand-structures`` action. This can be considered a command
line interface to :class:`ost.mol.alg.ligand_scoring.LigandScorer` and more
information about arguments and outputs can be found there.

Details on the usage (output of ``ost compare-ligand-structures --help``):

.. code-block:: console

  usage: ost compare-ligand-structures [-h] -m MODEL [-ml [MODEL_LIGANDS ...]]
                                       -r REFERENCE
                                       [-rl [REFERENCE_LIGANDS ...]] [-o OUTPUT]
                                       [-mf {pdb,cif,mmcif}]
                                       [-rf {pdb,cif,mmcif}] [-of {json,csv}]
                                       [-mb MODEL_BIOUNIT]
                                       [-rb REFERENCE_BIOUNIT] [-ft] [-rna]
                                       [-sm] [-cd COVERAGE_DELTA] [-v VERBOSITY]
                                       [--full-results] [--lddt-pli]
                                       [--lddt-pli-radius LDDT_PLI_RADIUS]
                                       [--lddt-pli-amc] [--rmsd]
                                       [--radius RADIUS]
                                       [--lddt-lp-radius LDDT_LP_RADIUS] [-fbs]
  
  Evaluate model with non-polymer/small molecule ligands against reference.
  
  Example: ost compare-ligand-structures \
      -m model.pdb \
      -ml ligand.sdf \
      -r reference.cif \
      --lddt-pli --rmsd
  
  Structures of polymer entities (proteins and nucleotides) can be given in PDB
  or mmCIF format.
  
  Ligands can be given as path to SDF files containing the ligand for both model
  (--model-ligands/-ml) and reference (--reference-ligands/-rl). If omitted,
  ligands will be detected in the model and reference structures. For structures
  given in mmCIF format, this is based on the annotation as "non polymer entity"
  (i.e. ligands in the _pdbx_entity_nonpoly mmCIF category) and works reliably.
  For structures given in legacy PDB format, this is based on the HET records
  which is usually only set properly on files downloaded from the PDB (and even
  then, this is not always the case). This is normally not what you want. You
  should always give ligands as SDF for structures in legacy PDB format.
  
  Polymer/oligomeric ligands (saccharides, peptides, nucleotides) are not
  supported.
  
  Only minimal cleanup steps are performed (remove hydrogens and deuteriums,
  and for structures of polymers only, remove unknown atoms and cleanup element
  column).
  
  Ligands in mmCIF and PDB files must comply with the PDB component dictionary
  definition, and have properly named residues and atoms, in order for
  ligand connectivity to be loaded correctly. Ligands loaded from SDF files
  are exempt from this restriction, meaning any arbitrary ligand can be assessed.
  
  Output can be written in two format: JSON (default) or CSV, controlled by the
  --output-format/-of argument.

  Without additional options, the JSON ouput is a dictionary with three keys:
  
   * "model_ligands": A list of ligands in the model. If ligands were provided
     explicitly with --model-ligands, elements of the list will be the paths to
     the ligand SDF file(s). Otherwise, they will be the chain name, residue
     number and insertion code of the ligand, separated by a dot.
   * "reference_ligands": A list of ligands in the reference. If ligands were
     provided explicitly with --reference-ligands, elements of the list will be
     the paths to the ligand SDF file(s). Otherwise, they will be the chain name,
     residue number and insertion code of the ligand, separated by a dot.
   * "status": SUCCESS if everything ran through. In case of failure, the only
     content of the JSON output will be "status" set to FAILURE and an
     additional key: "traceback".
  
  Each score is opt-in and the respective results are available in three keys:
  
   * "assigned_scores": A list with data for each pair of assigned ligands.
     Data is yet another dict containing score specific information for that
     ligand pair. The following keys are there in any case:
  
      * "model_ligand": The model ligand
      * "reference_ligand": The target ligand to which model ligand is assigned to
      * "score": The score
      * "coverage": Fraction of model ligand atoms which are covered by target
        ligand. Will only deviate from 1.0 if --substructure-match is enabled.
  
   * "model_ligand_unassigned_reason": Dictionary with unassigned model ligands
     as key and an educated guess why this happened.
  
   * "reference_ligand_unassigned_reason": Dictionary with unassigned target ligands
     as key and an educated guess why this happened.
  
  If --full-results is enabled, another element with key "full_results" is added.
  This is a list of data items for each pair of model/reference ligands. The data
  items follow the same structure as in "assigned_scores". If no score for a
  specific pair of ligands could be computed, "score" and "coverage" are set to
  null and a key "reason" is added giving an educated guess why this happened.

  CSV output is a table of comma-separated values, with one line for each
  reference ligand. The following column is always available:

   * reference_ligand: If reference ligands were provided explicitly with
     --reference-ligands, elements of the list will be the paths to the ligand
     SDF file(s). Otherwise, they will be the chain name, residue number and
     insertion code of the ligand, separated by a dot.

  If lDDT-PLI was enabled with --lddt-pli, the following columns are added:

   * "lddt_pli", "lddt_pli_coverage" and "lddt_pli_model_ligand" are the
     lDDT-PLI score result, the corresponding coverage and assigned model ligand,
     if an assignment was found, respectively, empty otherwise.
   * "lddt_pli_unassigned" is empty if an assignment was found, otherwise it
     lists the short reason this reference ligand was unassigned.

  If BiSyRMSD was enabled with --rmsd, the following columns are added:

   * "rmsd", "rmsd_coverage". "rmsd_lddt_lp" "rmsd_bb_rmsd" and
     "rmsd_model_ligand" are the BiSyRMSD, the corresponding coverage,
     lDDT-LP, backbone RMSD and assigned model ligand, if an assignment was
     found, respectively, empty otherwise.
   * "rmsd_unassigned" is empty if an assignment was found, otherwise it
     lists the short reason this reference ligand was unassigned.
  
  options:
    -h, --help            show this help message and exit
    -m MODEL, --mdl MODEL, --model MODEL
                          Path to model file.
    -ml [MODEL_LIGANDS ...], --mdl-ligands [MODEL_LIGANDS ...], --model-ligands [MODEL_LIGANDS ...]
                          Path to model ligand files.
    -r REFERENCE, --ref REFERENCE, --reference REFERENCE
                          Path to reference file.
    -rl [REFERENCE_LIGANDS ...], --ref-ligands [REFERENCE_LIGANDS ...], --reference-ligands [REFERENCE_LIGANDS ...]
                          Path to reference ligand files.
    -o OUTPUT, --out OUTPUT, --output OUTPUT
                          Output file name. Default depends on format: out.json or
                          out.csv
    -mf {pdb,cif,mmcif}, --mdl-format {pdb,cif,mmcif}, --model-format {pdb,cif,mmcif}
                          Format of model file. pdb reads pdb but also pdb.gz,
                          same applies to cif/mmcif. Inferred from filepath if
                          not given.
    -rf {pdb,cif,mmcif}, --reference-format {pdb,cif,mmcif}, --ref-format {pdb,cif,mmcif}
                          Format of reference file. pdb reads pdb but also
                          pdb.gz, same applies to cif/mmcif. Inferred from
                          filepath if not given.
    -of {json,csv}, --out-format {json,csv}, --output-format {json,csv}
                          Output format, JSON or CSV, in lowercase. default: json
    -mb MODEL_BIOUNIT, --model-biounit MODEL_BIOUNIT
                          Only has an effect if model is in mmcif format. By
                          default, the asymmetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the ID (as a string) of the one which
                          should be used.
    -rb REFERENCE_BIOUNIT, --reference-biounit REFERENCE_BIOUNIT
                          Only has an effect if reference is in mmcif format. By
                          default, the asymmetric unit (AU) is used for scoring.
                          If there are biounits defined in the mmcif file, you
                          can specify the ID (as a string) of the one which
                          should be used.
    -ft, --fault-tolerant
                          Fault tolerant parsing.
    -rna, --residue-number-alignment
                          Make alignment based on residue number instead of
                          using a global BLOSUM62-based alignment (NUC44 for
                          nucleotides).
    -sm, --substructure-match
                          Allow incomplete (ie partially resolved) target
                          ligands.
    -cd COVERAGE_DELTA, --coverage-delta COVERAGE_DELTA
                          Coverage delta for partial ligand assignment.
    -v VERBOSITY, --verbosity VERBOSITY
                          Set verbosity level. Defaults to 3 (INFO).
    --full-results        Outputs scoring results for all model/reference ligand
                          pairs and store as key "full_results"
    --lddt-pli            Compute lDDT-PLI scores and store as key "lddt_pli".
    --lddt-pli-radius LDDT_PLI_RADIUS
                          lDDT inclusion radius for lDDT-PLI.
    --lddt-pli-amc        Add model contacts (amc) when computing lDDT-PLI.
    --rmsd                Compute RMSD scores and store as key "rmsd".
    --radius RADIUS       Inclusion radius to extract reference binding site
                          that is used for RMSD computation. Any residue with
                          atoms within this distance of the ligand will be
                          included in the binding site.
    --lddt-lp-radius LDDT_LP_RADIUS
                          lDDT inclusion radius for lDDT-LP.
    -fbs, --full-bs-search
                          Enumerate all potential binding sites in the model
                          when searching rigid superposition for RMSD
                          computation
