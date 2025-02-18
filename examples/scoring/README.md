Structure comparison examples
=============================

OpenStructure provides "actions" for general benchmarking use cases.

* **compare-structures**: Focuses on comparisons involving polymer entities, i.e.
  protein, DNA and RNA chains
* **compare-ligand-structures**: Focuses on comparisons of interactions between
  polymer entities and non-polymer entities, i.e. small molecule ligands
  
The example commands here assume an OpenStructure installation
(compile instructions: https://openstructure.org/docs/install/). 
Running the computations in containers provide a considerably easier setup than
compiling OpenStructure from source. Instructions for setup and running
equivalent computations are available for

* [Docker](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/docker)
* [Singularity](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/singularity)

A detailed list of options can be found in the
[action documentation](https://openstructure.org/docs/actions/).
Alternatively, the command

```
ost <ACTION> -h
```

will list all available options of the respective action. Both actions compute
scores on an “opt-in” basis and produce output in JSON format. Example to compute
global and per-residue LDDT scores, as well as as QS-scores which are written to
default output (out.json):

```
ost compare-structures -m model.pdb -r reference.cif.gz --lddt --local-lddt --qs-score
```

An example output can be found [here](compare-structures_example_out.json) and
we refer to the action documentation for in-depth description of the provided
data items.

By default, model-reference chains are aligned using Needleman-Wunsch.
Many benchmarking efforts such as CASP and CAMEO assume residue numbers
according to target sequence(s). Both "actions" allow to derive model-reference
chain alignments according to these numbers which should be preferred in these
cases. This can be enabled by adding a `-rna` (residue number alignment) flag:

```
ost compare-structures -m model.pdb -r reference.cif.gz --lddt --local-lddt --qs-score -rna
```

The same example also contains small molecule ligands.
We can compute LDDT-PLI and BiSyRMSD with:

```
ost compare-ligand-structures -m model.pdb -r reference.cif.gz -ml *.sdf --rmsd --lddt-pli
```

An example output can be found [here](compare-ligand-structures_example_out.json)
and we refer to the action documentation for in-depth description of the provided
data items.

Again, it is advised to use the `-rna` flag if applicable. In this example,
reference ligands are directly extracted from the provided mmCIF file based on
"non-polymer" entity types.
This only works in case of mmCIF input AND if the respective ligand is in the
PDB component dictionary which defines connectivity (matching based on compound
name).
Container solutions come with such a dictionary which has been created at build
time. Check the Docker/Singularity instructions linked above on how to set the
latest dictionary, i.e. "Compound Library".
You can override automatic extraction by providing SDF files with ligand
coordinates and connectivity information. If the receptor is provided in
PDB format, ligands must be provided in SDF format.

