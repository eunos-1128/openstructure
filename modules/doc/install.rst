Installing OpenStructure From Source
================================================================================

Brief Overview
--------------------------------------------------------------------------------

For a simple and portable way to use OpenStructure we recommend using conda or a
container solution.  We provide a conda package and recipes to build images for
`Docker <https://www.docker.com/>`_ and
`Singularity <https://sylabs.io/singularity/>`_.
The latest recipes and instructions can be found on our
`GitLab site <https://git.scicore.unibas.ch/schwede/openstructure/>`_, including
a link to OpenStructure's own `GitLab Docker registry <https://git.scicore.unibas.ch/schwede/openstructure/container_registry>`_ (`Docker instructions`_ and
`Singularity instructions`_).

If you wish to compile OpenStructure outside of a container, you need to follow
the steps which we describe in detail below. In essence, these steps are:

* `Installing the Dependencies`_
* `Getting the Source Code`_
* `Configuring the build`_
* `Building the Project`_
* `Building the Compound Library`_
* `Running the tests`_
* `Installing OpenStructure`_


Installing the Dependencies
--------------------------------------------------------------------------------

OpenStructure requires a C++11 enabled compiler (e.g. recent gcc/clang) and uses 
a bunch of open-source libraries. If you haven't already installed them, please 
install them now! Where appropriate, the preferred version is given in 
parentheses.

* `CMake <http://cmake.org>`_ (3.25)
* `Python3 <http://python.org>`_ (3.9.5)
* `Boost <http://boost.org>`_ (1.76.0)
* `zlib <https://zlib.net/>`_ (usually comes with Boost or system)
* `Eigen3 <http://eigen.tuxfamily.org>`_ (3.4.0)
* `SQLite3 <https://www3.sqlite.org>`_ (3.35.4)
* `FFTW3 <http://fftw.org>`_ (3.3.9). By default, OpenStructure is compiled with single
  precision and thus also requires FFTW to be compiled with single precision.
  Most platforms offer this as a second package. If you are compiling manually,
  use the `--enable-single` option.
* `libtiff <http://www.libtiff.org>`_ (4.2.0)
* `libpng <http://www.libpng.org>`_ (1.6.37, also needed for GUI)

If you would like to use the info module, also install:

* `Qt5 <http://qt-project.org/>`_ 

If you would like to use the graphical user interface (GUI), also install:

* `Qt5 <http://qt-project.org/>`_ 
* `SIP <http://www.riverbankcomputing.co.uk/software/sip/download>`_
* `PyQt5 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`_

If you would like to use the :mod:`molecular mechanics <ost.mol.mm>` module:

* `OpenMM <https://simtk.org/home/openmm>`_ (7.7.0)

If you would like to do pairwise sequence alignments with parasail
in the :mod:`bindings <ost.bindings>` module:

* `parasail <https://github.com/jeffdaily/parasail/>`_ (2.6.2)

We do not provide backwards compatibility to Python 2.7. The last
release supporting Python 2.7 is 1.11.0.


Getting the Source Code
--------------------------------------------------------------------------------

OpenStructure uses `git` as the revision control system. The main repository can
be browsed `here <https://git.scicore.unibas.ch/schwede/openstructure.git>`_. To
get the source code, use git clone:

.. code-block:: bash

  git clone https://git.scicore.unibas.ch/schwede/openstructure.git <directory-name>
  
The above command will clone OpenStructure into the directory specified by
`<directory-name>`. If omitted, the directory will be called openstructure. 

.. note::

  Some versions of curl have trouble with the certificate of the OpenStructure
  git server and fail to clone the repository. To work around this, disable the
  SSL certificate verification by setting the following environment variable:
  
  .. code-block:: bash

    export GIT_SSL_NO_VERIFY=1


Picking the right branch
--------------------------------------------------------------------------------

By default you are checking out the master branch. Master is by definition a
stable branch. It always points to the latest release. However, there are
several other branches at your disposal. The main development is happening in
the develop branch. It contains the newest features and bug fixes. However, we
don't make any guarantees that the develop branch is bug free and doesn't
contain major bugs. After all, it's in constant flux. If you are developing new
features, start your feature branch off develop. Besides that, there are several
smaller features branches that are used to group together commits for one
specific features. To change to a specific branch, use

.. code-block:: bash

  git checkout <branch-name>


Configuring the build
--------------------------------------------------------------------------------

OpenStructure uses `CMake <http://cmake.org>`_ for compiling and building the
project. The next required step is to configure the build environment using
cmake. You can do that by invoking `cmake` in the project directory.

.. code-block:: bash

  cmake . <options>

There are two kinds of options: Options that let you control the building
behaviour, enabling and disabling the compilation of certain modules and options
that let you tell CMake where to find the dependencies. All of them are passed
to CMake via `-D<opt>=<value>`.


.. _cmake-flags:

Flags to Control the Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, `CMake <http://cmake.org>`_ searches the standard directories for
dependencies. However, on some systems, this might not be enough. Here is a
short description of how CMake figures out what dependencies to take and how you
can influence it.

* Boost is mainly controlled via the `BOOST_ROOT` option. If boost wasn't
  found, it should be set to the prefix of the boost installation. If for some
  reason, it is desirable to use the non-multithreaded boost libraries, you can
  switch `Boost_USE_MULTITHREADED` off (it is on by default).

* `Python_ROOT_DIR` is the Python equivalent of BOOST_ROOT. It should be set to 
  the prefix path containing the python binary, headers and libraries.

* `SYS_ROOT` controls the general prefix for searching libraries and headers.
  By default, it is set to `/`.
  
* `COMPOUND_LIB` specifies the location of the compound library and
  activates the rule-based-builder. The compound library is based on 
  the component dictionary released by the PDB, and it specifies atoms
  of a certain residue or connectivities between atoms etc. The 
  :doc:`compound library <conop/compoundlib>` itself is created from the 
  component dictionary by calling the OpenStructure chemdict_tool. 
  By default this is switched off but it is highly recommended to provide a
  compound library to use all features of OpenStructure.

* `ENABLE_GUI` controls whether to build the graphical user interface module.
  By default, this is switched on.

* `ENABLE_GFX` controls whether to build the graphics module. By default, this
  is switched on. If it is switched off, it also switches `ENABLE_GUI` off.

* `ENABLE_INFO` controls whether to build the info module. By default, this is
  switched on. If it is switched off, it also switches `ENABLE_GFX` off and
  removes all dependencies to Qt.

* `QT_QMAKE_EXECUTABLE` defines the exact Qt installation to take. It should 
  be set to the full path to `qmake`. This is only needed if `ENABLE_INFO` is
  switched on.

* `COMPILE_TMTOOLS` will activate bindings for TMAlign and TMScore, which are 
  then available at python level. This option requires a Fortran compiler. 
  By default, this option is switched off.

* `ENABLE_MM` controls whether the molecular mechanics module is enabled. By
  default, this is switched off. If it is turned on, you should also set the
  paths to your local OpenMM installation:

  * `OPEN_MM_INCLUDE_DIR`: the include path
  * `OPEN_MM_LIBRARY`: the libOpenMM library
  * `OPEN_MM_PLUGIN_DIR`: the path for OpenMM plugins
  * see example below for commonly used paths

* `ENABLE_PARASAIL` controls whether parasail should be enabled. By default,
  this is switched off. If it is turned on, you must also set the paths
  to your parasail installation:

  * `PARASAIL_INCLUDE_DIR`: the include path containing the file parasail.h
  * `PARASAIL_LIBRARY`: the parasail library

* Several paths to other libraries can be set if they are not in the expected
  locations:

  * `Python_LIBRARY` defines the location of the Python library (file name
    starting with `libpython`). This must be set if it is not in
    `$Python_ROOT_DIR/lib`.
  * `EIGEN3_INCLUDE_DIR` defines the include folder of Eigen3 (contains `Eigen`
    folder with include files).
  * `FFTW_LIBRARY` defines the location of the FFTW3 library (file name starting
    with `libfftw3f` (or `libfftw3` if `USE_DOUBLE_PRECISION` is switched on))
  * `FFTW_INCLUDE_DIR` defines the include folder of FFTW3 (contains include
    files directly)
  * `PNG_LIBRARY` defines the location of the libpng library (file name starting
    with `libpng`)
  * `PNG_PNG_INCLUDE_DIR` defines the include folder of libpng (contains include
    files directly)
  * `ZLIB_LIBRARY` defines the location of the zlib library (file name starting
    with `libz`)
  * `ZLIB_INCLUDE_DIR` defines the include folder of zlib (contains include
    files directly)
  * `TIFF_LIBRARY` defines the location of the libtiff library (file name
    starting with `libtiff`)
  * `TIFF_INCLUDE_DIR` defines the include folder of libtiff (contains include
    files directly)
  * `SQLITE3_LIBRARY` defines the location of the SQLite3 library (file name starting
    with `libsqlite3`)
  * `SQLITE3_INCLUDE_DIR` defines the include folder of SQLite3 (contains include
    files directly)
  * Usually, you will receive errors for those variables when executing `cmake`
    and set them accordingly as needed.

* `OPENGLPREFERENCE_LEGACY` switches the GL implementation to be used by OpenGL.
  The default is what should be used on modern systems. But since there are some
  reports on the internet claiming that the default does not work everywhere,
  this switch enables the usage of the legacy implementation of GL.
  
Build Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `OPTIMIZE` can be switched on to build an optimised (-O3 -DNDEBUG) version of
  OpenStructure. By default, this is switched off.

* `USE_DOUBLE_PRECISION` will switch on double precision within OpenStructure. 
  By default, this is switched off.

* `ENABLE_STATIC` allows some parts of OpenStructure to be statically linked 
  and thus can be used more easily across a heterogeneous setup, e.g. older 
  systems and newer systems. Note that enabling this flag will not compile the
  full OpenStructure package and it is not guaranteed to lead to fully portable
  binaries. By default, this is switched off.

* For deployment of OpenStructure with `make install` there are two relevant
  settings to consider:

  * `PREFIX` or `CMAKE_INSTALL_PREFIX` are used to define the path where the
    OpenStructure `stage` folder will be installed to.
  * `USE_RPATH` can be switched on to embed rpath upon make install. By default,
    this option is switched off.

* Experimental settings (only change if you know what you are doing):

  * `USE_SHADER` controls whether to compile with shader support. By default,
    this is turned off.
  * `ENABLE_SPNAV` controls whether 3DConnexion devices should be supported. By
    default, this is turned off.
  * `PROFILE` can be switched on to enable a (very verbose) code profiler. By
    default, this is turned off.
  * `UBUNTU_LAYOUT` can be turned on to switch the directory layout of the
    `stage` folder to be more ubuntu-like. By default, this is switched off.
  * `HIDDEN_VISIBILITY` can be turned on to add "-fvisibility=hidden" to gcc's
    compile flags (only if GNU compiler used). By default, this is switched off.

Known Issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Depending on how the dependecies (e.g. Boost) are compiled, linking might fail
  with something like: `error: undefined reference to pthread_condattr_destroy`.
  Add "-pthread" to the linking options by appending the following to your cmake
  command: `-DCMAKE_EXE_LINKER_FLAGS=" -pthread"`

Example Configurations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Generic Linux without GUI**

The simplest way to compile OpenStructure is to disable the GUI and any
dependency to Qt5. You can build an optimised OpenStructure without GUI as
follows:

.. code-block:: bash

  cmake . -DOPTIMIZE=ON -DENABLE_INFO=OFF

The molecular mechanics module can be enabled by installing OpenMM and adding
the appropriate flags as follows (replace `<OPENMM>` with the actual path to
OpenMM):

.. code-block:: bash

  cmake . -DOPTIMIZE=ON -DENABLE_INFO=OFF -DENABLE_MM=ON \
          -DOPEN_MM_LIBRARY=<OPENMM>/lib/libOpenMM.so \
          -DOPEN_MM_INCLUDE_DIR=<OPENMM>/include/ \
          -DOPEN_MM_PLUGIN_DIR=<OPENMM>/lib/plugins

Note that the OpenMM binaries available online may be incompatible with files
compiled using your gcc compiler (known as "Dual ABI" issue). This has been
observed for OpenMM versions 6.1 until 7.1.1 when compiling with gcc versions >=
5.1. In those cases, you cannot use the binaries and will have to install OpenMM
from source.


**Ubuntu 24.04 LTS**

Besides the molecular mechanics module, we also enable parasail here.
All the dependencies can be installed from the package manager as follows:

.. code-block:: bash

  sudo apt-get install cmake g++ libtiff-dev libfftw3-dev libeigen3-dev \
               libpng-dev python3-all python3-pyqt5 libboost-all-dev \
               qt5-qmake qtbase5-dev libpng-dev libsqlite3-dev \
               libopenmm-dev libopenmm-plugins libparasail-dev

Now, all dependencies are located in standard locations and cmake will
automatically find them. As a single quirk, we need to specify the
OpenMM plugin directory. Lets do a proper out of source build here:

.. code-block:: bash

  mkdir build
  cd build
  cmake .. -DOPTIMIZE=ON -DENABLE_MM=1 -DENABLE_PARASAIL=1 \
           -DOPEN_MM_PLUGIN_DIR=/lib/x86_64-linux-gnu/openmm/plugins

Building the Project
--------------------------------------------------------------------------------

Type ``make``. If you are using a multi-core machine, you can use the `-j` flag
to run multiple jobs at once.

.. code-block:: bash

  make

Building the Compound Library
--------------------------------------------------------------------------------

One thing is missing for a fully functional OpenStructure installation.
The compound library. It is used at various places for connectivity
information and certain algorithms do not work without.
Besides an OpenStructure executable, we just built the
:ref:`chemdict_tool <mmcif-convert>` which converts the PDB chemical component dictionary
into our internal format:

.. code-block:: bash

  wget https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
  stage/bin/chemdict_tool create components.cif.gz compounds.chemlib -i

We can rerun cmake and make. All cmake parameters from the original
configuration remain in the cache.

.. code-block:: bash

  cmake .. -DCOMPOUND_LIB=compounds.chemlib
  make

Running the tests
--------------------------------------------------------------------------------

Many parts of OpenStructure are covered by unit tests. You can run them with:

.. code-block:: bash

  make check


Building the documentation
--------------------------------------------------------------------------------

Many parts of OpenStructure are covered by unit tests. You can run them with:

.. code-block:: bash

  stage/bin/ost doc/make.py


Installing OpenStructure
--------------------------------------------------------------------------------

To make OpenStructure available system-wide, the easiest is to install it with:

.. code-block:: bash

  make install

This will install the OpenStructure binaries and libraries to the prefix
folder specified by `CMAKE_INSTALL_PREFIX` in the cmake configuration step
(usually this will be `/usr/local`).

In many distributions (such as Ubuntu), libraries are not picked up 
automatically from `/usr/local` so you will need to export a few environment 
variables to make things work:

.. code-block:: bash

  export OST_ROOT="/usr/local"
  export PYTHONPATH="/usr/local/lib64/python3.10/site-packages"
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib64:/usr/local/openmm/lib/"


If you don't want to install OpenStructure system-wide, for instance if this is
a development build, you can also configure environment variables to point 
directly to the `stage` folder.

.. code-block:: bash

  export OST_ROOT=/path/to/ost/stage
  export PATH=$OST_ROOT/bin:$PATH
  export PYTHONPATH=$OST_ROOT/lib64/python3.10/site-packages/:$PYTHONPATH


What's next?
--------------------------------------------------------------------------------

On Linux and macOS, you can start dng from the command-line. The binaries are
all located in stage/bin:

.. code-block:: bash

  dng
  
or, to start the command-line interpreter:

.. code-block:: bash

  ost

The `ost` executable can run python scripts. For instance we can run the Docker
test script:

.. code-block:: bash

  ost docker/test_docker.py

This should produce some output and announce that "OST is working!"

You can also import OpenStructure directly into your existing python scripts,
jupyter notebooks etc. Simply make sure to point the following environment
variables to the right folders:

.. code-block:: bash

  export OST_ROOT=/path/to/ost/stage
  export PYTHONPATH=$OST_ROOT/lib64/python3.10/site-packages/:$PYTHONPATH
  python

And then you can simply import ost as a module:

.. code-block:: python

  import ost


..  LocalWords:  Homebrew cmake CMake zlib SQLite FFTW libtiff libpng PyQt OST
..  LocalWords:  SSL macOS Makefiles PDB qmake PNG libz libsqlite OPTIMIZE dng
..  LocalWords:  DNDEBUG RPATH rpath SHADER shader SPNAV DConnexion profiler
..  LocalWords:  DOPTIMIZE DENABLE DOPEN DPYTHON DBOOST DSYS Xcode Eigen Sur
..  LocalWords:  Monterey SDKROOT DPython DIRS DCMAKE isystem CXX
