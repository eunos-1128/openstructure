# OpenStructure - A computational structural biology framework

OpenStructure provides a modular, flexible, molecular modelling environment
which allows to integrate, process and visualize information of different origin
such as sequences, alignments and 3D structures.

Please refer to www.openstructure.org for more information and documentation.

Thank you for you interest and enjoy the straightforward way of handling protein
structure data!

Please do not hesitate to contact us for feedback or troubleshooting:

 openstructure-users@maillist.unibas.ch


## OpenStructure Installation

OpenStructure is developed and tested across various Linux distributions. You can find detailed build instructions and a list of required dependencies here:
https://openstructure.org/docs/dev/install/

### Prefer Not to Compile it Yourself?

For a simple and portable setup, we recommend using a containerized
solution. OpenStructure provides its own Docker container registry,
making deployment easier.

### Containerized Solutions:

* Docker: [OpenStructure Docker Instructions](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/docker)
* Singularity: [OpenStructure Singularity Instructions](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/singularity)

## Getting started

The installation and container instructions will help you set up
OpenStructure for scripting. Additionally, here are some advanced
use cases worth highlighting.

### Benchmarking and Scoring

OpenStructure implements a benchmarking suite for comparing macromolecular
complexes. To get started, check out the examples in the `examples/scoring`
directory of this repository.

### Modeling

ProMod3 is a fully featured modeling engine built on top of OpenStructure.
It powers the SWISS-MODEL web server and provides advanced modeling
capabilities. Learn more here:

[ProMod3 Documentation](https://openstructure.org/promod3)


## Cite Us

If you like our software and have used it in your research project, please cite
the following paper on OpenStructure:

 M. Biasini, T. Schmidt, S. Bienert, V. Mariani, G. Studer, J. Haas, N. Johner,
 A.D. Schenk, A. Philippsen and T. Schwede, OpenStructure: an integrated
 software framework for computational structural biology, Acta Cryst., 2013

If you use the code or binary in OpenStructure to compute lDDT scores, please
also cite the following reference:

 V. Mariani, M. Biasini, A. Barbato, T. Schwede, lDDT : A local superposition-
 free score for comparing protein structures and models using distance
 difference tests, Bioinformatics, 2013

If you use the code in OpenStructure to compute QS scores, please also cite the
following reference:

 M. Bertoni, F. Kiefer, M. Biasini, L. Bordoli, T. Schwede, Modeling protein
 quaternary structure of homo- and hetero-oligomers beyond binary interactions
 by homology, Scientific Reports, 2017 

================= The OpenStructure Team =======================================
