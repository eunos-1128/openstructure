# OpenStructure - A computational structural biology framework

OpenStructure provides a modular, flexible, molecular modelling environment
which allows to integrate, process and visualize information of different origin
such as sequences, alignments and 3D structures.

Please refer to www.openstructure.org for more information and documentation.

Thank you for your interest and enjoy the straightforward way of handling protein
structure data!

Please do not hesitate to contact us for feedback or troubleshooting:

 <a href="mailto:&#111;&#112;&#101;&#110;&#115;&#116;&#114;&#117;&#099;&#116;&#117;&#114;&#101;&#045;&#117;&#115;&#101;&#114;&#115;&#064;&#109;&#097;&#105;&#108;&#108;&#105;&#115;&#116;&#046;&#117;&#110;&#105;&#098;&#097;&#115;&#046;&#099;&#104;">&#111;&#112;&#101;&#110;&#115;&#116;&#114;&#117;&#099;&#116;&#117;&#114;&#101;&#045;&#117;&#115;&#101;&#114;&#115;&#064;&#109;&#097;&#105;&#108;&#108;&#105;&#115;&#116;&#046;&#117;&#110;&#105;&#098;&#097;&#115;&#046;&#099;&#104;</a>

## OpenStructure Installation

For a simple and portable setup, we recommend using a containerized
solution. OpenStructure provides its own Docker container registry,
making deployment easier. Deploying a docker image just needs a
docker pull which typically takes less than a minute. Singularity
containers bootstrap from the docker container but must be built
by the user. Detailed instructions can be found here:

* Docker: [OpenStructure Docker Instructions](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/docker)
* Singularity: [OpenStructure Singularity Instructions](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/singularity)

OpenStructure is developed and tested across various Linux distributions.
You can find detailed build instructions and a list of required dependencies here:
https://openstructure.org/docs/install/

## Getting started

The installation and container instructions will help you set up
OpenStructure for scripting. Additionally, here are some advanced
use cases worth highlighting.

### Benchmarking and Scoring

OpenStructure implements a benchmarking suite for comparing macromolecular
complexes. To get started, check out the examples in the [examples/scoring](examples/scoring)
directory of this repository.

### Modeling

[ProMod3](https://git.scicore.unibas.ch/schwede/ProMod3) is a fully featured modeling engine built on top of OpenStructure.
It powers the [SWISS-MODEL](https://swissmodel.expasy.org) web server and provides advanced modeling
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

