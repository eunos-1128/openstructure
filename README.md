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

### Conda

OpenStructure is available as a [conda package on
Bioconda](https://bioconda.github.io/recipes/openstructure/README.html) for a simple and portable setup.
You can install it with the following command:

```
conda install bioconda::openstructure
```

Tested with [miniforge](https://conda-forge.org/miniforge/).

You can then run a test script:

```
wget https://git.scicore.unibas.ch/schwede/openstructure/-/raw/master/docker/test_docker.py -O test_ost.py
ost test_ost.py
```

### Containers

For a fully portable setup, we provide containerized solutions.
OpenStructure provides its own Docker container registry,
making deployment easier. Deploying a docker image just needs a
docker pull which typically finishes in about a minute depending
on your local hardware and internet connection. Singularity
containers bootstrap from the docker container but must be built
by the user. Both solutions require root permissions.


#### Docker

For complete documentation on using Docker with OpenStructure, 
[click here](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/docker).
This section provides a quick-start guide to help you get OpenStructure up and running using Docker.
Most docker installations require you to add `sudo` in front of the docker commands.

Get the latest Docker image from the OpenStructure registry:

```
docker pull registry.scicore.unibas.ch/schwede/openstructure:latest
```

And run a test script:

```
wget https://git.scicore.unibas.ch/schwede/openstructure/-/raw/master/docker/test_docker.py -O test_ost.py
docker run --rm -v $(pwd):/home registry.scicore.unibas.ch/schwede/openstructure:latest test_ost.py
```

#### Singularity

For complete documentation on using Singularity with OpenStructure, 
[click here](https://git.scicore.unibas.ch/schwede/openstructure/tree/master/singularity).
This section provides a quick-start guide to help you get OpenStructure up and running using Singularity.

Building the singularity container requires root permissions:

```
wget https://git.scicore.unibas.ch/schwede/openstructure/-/raw/master/singularity/Singularity
sudo singularity build ost.img Singularity
```

And run a test script:

```
wget https://git.scicore.unibas.ch/schwede/openstructure/-/raw/master/docker/test_docker.py -O test_ost.py
singularity run --app OST ost.img test_ost.py
```

### Build from source

OpenStructure is developed and tested across various Linux distributions.
You can find [detailed build instructions and a list of required dependencies
in the online documentation](https://openstructure.org/docs/install/).

## Getting started

Once you have OpenStructure installed, check out the [online
documentation](https://openstructure.org/docs/).

Additionally, here are some advanced use cases worth highlighting.

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

