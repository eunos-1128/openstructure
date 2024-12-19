FROM ubuntu:22.04

# ARGUMENTS
###########
ARG OPENSTRUCTURE_VERSION="2.9.1"
ARG SRC_FOLDER="/usr/local/src"
ARG CPUS_FOR_MAKE=2

ARG DEBIAN_FRONTEND=noninteractive

# INSTALL SYSTEM DEPS
#####################
RUN apt-get update -y && apt-get install -y cmake \
                                            g++ \
                                            gfortran \
                                            wget \
                                            tar \
                                            libsqlite3-dev \
                                            sip-dev \
                                            libtiff-dev \
                                            libfftw3-dev \
                                            libeigen3-dev \
                                            libboost-all-dev \
                                            libpng-dev \
                                            python3-all \
                                            python3-numpy \
                                            python3-scipy \
                                            python3-pandas \
                                            python3-networkx \
                                            clustalw \
                                            voronota \
                                            libopenmm-dev \
                                            libopenmm-plugins \
                                            libparasail-dev \
                                            locales && \
                                            # CLEANUP
                                            rm -rf /var/lib/apt/lists/*
            
# INSTALL OST
#############
RUN cd ${SRC_FOLDER} && \
    # copy ost release
    wget -O openstructure-${OPENSTRUCTURE_VERSION}.tar.gz -nc https://git.scicore.unibas.ch/schwede/openstructure/-/archive/${OPENSTRUCTURE_VERSION}/openstructure-${OPENSTRUCTURE_VERSION}.tar.gz && \
    mkdir openstructure-${OPENSTRUCTURE_VERSION} && \
    tar xf openstructure-${OPENSTRUCTURE_VERSION}.tar.gz -C ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION} --strip-components=1 && \
    mkdir -p ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build && \
    cd ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build && \
    cmake .. -DOPTIMIZE=1 \
             -DENABLE_MM=1 \
             -DOPEN_MM_PLUGIN_DIR=/lib/x86_64-linux-gnu/openmm/plugins \
             -DENABLE_PARASAIL=1 \
             -DCOMPILE_TMTOOLS=1 \
             -DENABLE_GFX=1 \
             -DENABLE_GUI=0 \
             -DENABLE_INFO=0 && \
    make -j ${CPUS_FOR_MAKE} && \
    wget https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz && \
    stage/bin/chemdict_tool create components.cif.gz compounds.chemlib pdb -i && stage/bin/chemdict_tool update ../modules/conop/data/charmm.cif compounds.chemlib charmm && \
    cmake .. -DCOMPOUND_LIB=${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}/build/compounds.chemlib && \
             make -j ${CPUS_FOR_MAKE} && make check && make install && \
    rm ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}.tar.gz && \
    rm -rf ${SRC_FOLDER}/openstructure-${OPENSTRUCTURE_VERSION}

# ENVIRONMENT
#############
WORKDIR /home
ENV OST_ROOT="/usr/local"
ENV PYTHONPATH="/usr/local/lib64/python3.10/site-packages"
ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib64:/usr/local/openmm/lib/"
ENV OPENSTRUCTURE_VERSION=$OPENSTRUCTURE_VERSION

ENTRYPOINT ["ost"]
CMD ["-i"]

