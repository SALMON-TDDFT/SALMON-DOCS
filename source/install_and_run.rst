.. _install-and-run:

Install and Run
================

Prerequisites
----------------

In this guide, it is assumed that readers have a basic knowledge of Linux and its command line operations.
For the installation of SALMON, following packages are required.

- Fortran90/C compiler. SALMON assumes users have one of the following compilers:

  - GCC (Gnu Compiler Collection)
  - Intel Compiler
  - Fujitsu Compiler (at FX100 / K-Computer)

- One of the following library packages for linear algebra:

  - Netlib BLAS/LAPACK/ScaLAPACK
  - Intel Math Kernel Library (MKL)
  - Fujitsu Scientific Subroutine Library 2 (SSL-II)

- Build tools:

  - CMake

If you use other compilers, you may need to change build scripts (CMake). See :any:`additional-options-in-configure`.
If no numerical library is installed on your computer system, you may need to install BLAS/LAPACK by yourself.
See :any:`troubleshooting-install`.

For the installation of SALMON, we adopt the CMake tools as the first option.
If there were any problems to use CMake tools in your environment, you may use the GNU make tools.
See :any:`troubleshooting-install`.

Download
-----------------

The newest version of SALMON can be downloaded from `download page <http://salmon-tddft.jp/download.html>`__.
To extract files from the downloaded file ``SALMON-<VERSION>.tar.gz``, type the following command in the command-line::

  $ tar -zxvf ./salmon-<VERSION>.tar.gz

After the extraction, the following directories will be created::

  SALMON
    |- src          Source codes
    |- example      Samples
    |- cmakefiles   CMake related files
    |- gnumakefiles GNU Makefiles for building


Build and Install
------------------

To compile SALMON to create executable the binary files, we adopt to use CMake tools as the first option.
In case you fail to build SALMON using CMake in your environment, we may use Gnu Make. See :any:`build-gnu-make`.


Checking CMake availability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, examine whether CMake is usable in your environment or not.
Type the following in Linux command-line::

    $ cmake --version

If CMake is not installed in your system, an error message such as ``cmake: command not found`` will appear.
If CMake is installed on your system, the version number will be shown.
To build SALMON, CMake of version 3.14.0 or later is required.
If you confirm that CMake of version 3.14.0 or later is installed in your system, proceed to :any:`build-cmake`.
However, we realize that old versions of CMake are installed in many systems.
If CMake is not installed or CMake of older versions is installed in your system, you need to install the new version by yourself.
It is a simple procedure and explained below.


Installation of CMake (pre-compiled binary of Linux)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CMake <https://cmake.org/>`_ is a cross-platform build tool.
The simplest way to make CMake usable in your environment is to get `the binary distribution of CMake from the download page <https://cmake.org/download/>`_. (The file name of the binary distribution will be ``cmake-<VERSION>-<PLATFORM>.tar.gz``). In standard Linux environment, a file for the platform of Linux x86_64 will be appropriate.

To download the file, proceed as follows: We assume that you are in the directory that you extracted files from the downloaded file of SALMON,
and that you will use the version 3.16.8. First get the URL of the download link from your browser, and use ``wget`` command in your Linux command-line::

    $ wget https://cmake.org/files/v3.16/cmake-3.16.8-Linux-x86_64.tar.gz

Next, unpack the archive by::

    $ tar -zxvf cmake-3.16.8-Linux-x86_64.tar.gz

and you will have the binary ``make-3.16.8-Linux-x86_64/bin/cmake`` in your directory.

To make the ``cmake`` command usable in your command-line, you need to modify the environment variable ``$PATH`` so that the executable of CMake are settled inside the directory specified in your ``$PATH``.
If you use the bash shell, you need to modify the file ``~/.bashrc`` that specifies the ``$PATH`` variable. It can be done by typing the following command in your login directory::

    $ export PATH=<SALMON_INSTALLATION_DIRECTORY>/cmake-3.16.8-Linux-x86_64/bin:$PATH

and then reload the configuration by typing::

    $ source ~/.bashrc

See :any:`installation-cmake` describes Other way of the installation.


.. _build-cmake:

Build using CMake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Confirming that CMake of version 3.14.0 or later can be usable in your environment, proceed the following steps.
We assume that you are in the directory SALMON.

1. Create a new temporary directory ``build`` and move to the directory::

    $ mkdir build
    $ cd build


2. Execute the python script ''configure.py'' and then make::

    $ python ../configure.py --arch=ARCHITECTURE --prefix=../
    $ make
    $ make install


In executing the python script, you need to specify ``ARCHITECTURE`` that indicates the architecture of the CPU in your computer system such as ``intel-avx``. The options of the ``ARCHITECUTRE`` are as follows:

================  =======================================  ================  =================
arch              Detail                                   Compiler          Numerical Library
================  =======================================  ================  =================
intel-knl         Intel Knights Landing                    Intel Compiler    Intel MKL
intel-knc         Intel Knights Corner                     Intel Compiler    Intel MKL
intel-avx         Intel Processer (Ivy-, Sandy-Bridge)     Intel Compiler    Intel MKL
intel-avx2        Intel Processer (Haswell, Broadwell ..)  Intel Compiler    Intel MKL
intel-avx512      Intel Processer (Skylake-SP)             Intel Compiler    Intel MKL
fujitsu-fx100     FX100 Supercomputer                      Fujitsu Compiler  SSL-II
fujitsu-a64fx-ea  A64FX processor (Fugaku, FX1000, FX700)  Fujitsu Compiler  SSL-II
================  =======================================  ================  =================

If the build is successful, you will get a file ``salmon`` at the top-level build directory.


Files necessary to run SALMON
------------------------------------

To run SALMON, at least two kinds of files are required for any calculations.
One is an input file with the filename extension ``*.inp`` that should be read from the standard input ``stdin``.
This file should be prepared in the Fortran90 namelist format.
Pseudopotential files of relevant elements are also required.
Depending on your purpose, some other files may also be necessary.
For example, coordinates of atomic positions of the target material may be either written in the input file or prepared as a separate file.


Pseudopotentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SALMON utilizes norm-conserving pseudpotentials.
You may find pseudopotentials of some elements in the samples prepared in :any:`Exercises`.
In SALMON, several formats of pseudopotentials may be usable.
Pseudopotentials with an extension ``.fhi`` can be obtained from the website listed below.
(This is a part of previous atomic data files for the ABINIT code.)

====================================  =====================================================================================
Pseudopotential                       Website
====================================  =====================================================================================
Pseudopotentials for the ABINIT code  https://www.abinit.org/sites/default/files/PrevAtomicData/psp-links/psp-links/lda_fhi
====================================  =====================================================================================

Filenames of the pseudopotentials should be written in the input file.


input file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Input files are composed of several blocks of namelists::

   &namelist1
     variable1 = int_value
     variable2 = 'char_value'
   /
   &namelist2
     variable1 = real8_value
     variable2 = int_value1, int_value2, int_value3
   /

A block of namelists starts with ``&namelist`` line and ends with ``/`` line.
The blocks may appear in any order.

Between two lines of ``&namelist`` and ``/``, descriptions of variables and their values appear.
Note that many variables have their default values so that it is not necessary to give values for all variables.
Descriptions of the variables may appear at any position if they are between ``&namelist`` and ``/``.

SALMON describes electron dynamics in systems with both isolated and periodic boundary conditions.
The boundary condition is specified by the variable ``iperiodic`` in the namelist ``&system``.

Calculations are usually achieved in two steps; first, the ground state calculation is carried out and then electron dynamics calculations in real time is carried out. A choice of the calculation mode or theory in the calculation is specified by the variable ``theory`` in the namelist ``&calculation``.
In the typical way, the ground state calculation based on DFT is first carried out specifying ``theory = 'dft'``.
Then the real-time electron dynamics calculation based on TDDFT is carried out specifying ``theory = 'tddft_pulse'``.

In :any:`Exercises`, we prepare six exercises that cover typical calculations feasible by SALMON.
We prepare explanations of the input files of the exercises that will help to prepare input files of your own interests.

There are more than 20 groups of namelists. A complete list of namelist variables is given in the file ``SALMON/manual/input_variables.md``.
Namelist variables that are used in our exercises are explained at :any:`Inputs`.


Run SALMON
-----------------------------------

Before running SALMON, the following preparations are required as described above: The executable file of ``salmon`` should be built from the source file of SALMON. An input file ``inputfile.inp`` and pseudopotential files should also be prepared.

The execution of the calculation can be done as follows: In single process environment, type the following command::

    $ salmon < inputfile.inp > fileout.out

In multiprocess environment in which the command to execute parallel calculations using MPI is ``mpiexec``, type the following command::

    $ mpiexec -n NPROC salmon < inputfile.inp > fileout.out

where NPROC is the number of MPI processes that you will use.

The execution command and the job submission procedure depends much on local environment. We summarize general conditions to execute SALMON:

- SALMON runs in both single-process and multi-process environments using MPI.
- executable files are prepared as ``salmon`` in the standard build procedure.
- to start calculations, ``inputfile.inp`` should be read through ``stdin``.


Appendix
------------

.. _additional-options-in-configure:

Additional options in configure.py script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Manual specifications of compiler and environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In executing ``configure.py``, you may manually specify compiler and environment variables instead of specifying the architecture, for example::

    $ python ../configure.py FC=mpiifort CC=mpiicc FFLAGS="-xAVX" CFLAGS="-restrict -xAVX"

The major options of ``configure.py`` are as follows:

=======================================  ===================================================
Commandline switch                       Detail
=======================================  ===================================================
-a ARCH, --arch=ARCH                     Target architecture
--enable-mpi, --disable-mpi              enable/disable MPI parallelization
--enable-scalapack, --disable-scalapack  enable/disable computations with ScaLAPACK library
--enable-eigenexa, --disable-eigenexa    enable/disable computations with RIKEN R-CCS EigenExa library
--enable-libxc, --disable-libxc          enable/disable computations with Libxc library
--with-lapack                            specified LAPACK/ScaLAPACK installed directory
--with-libxc                             specified Libxc installed directory
--debug                                  enable debug build
--release                                enable release build
FC, FFLAGS                               User-defined Fortran Compiler, and the compiler options
CC, CFLAGS                               User-defined C Compiler, and the compiler options
=======================================  ===================================================

In the build procedure by CMake, they search the following libraries.
If the libraries don't found in the path that is specified by environment variables, they will build the required libraries automatically.

- Netlib LAPACK (includes BLAS), and ScaLAPACK

    - We will download and build the Netlib libraries as the typical implementation.
    - http://www.netlib.org/lapack/
    - http://www.netlib.org/scalapack/

- Libxc

    - https://www.tddft.org/programs/libxc/

EigenExa will download and build automatically even if the library is installed to your machine.


Build for single process calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the ``--arch`` option, MPI parallelization is enabled as default in almost case.
If you use a single processor machine, explicitly specify ``--disable-mpi`` in executing the python script::

    $ python ../configure.py --arch=<ARCHITECTURE> --disable-mpi


Build by user-specified compiler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want that specify the compiler, set the ``FC`` and ``CC`` flags in executing the python script::

    $ python ../configure.py FC=gfortran CC=gcc

When not using the ``--arch`` option, MPI parallelization is disabled as default.



.. _build-gnu-make:

Build using GNU Makefile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If CMake build fails in your environment, we recommend you to try to use Gnu Make for the build process.
First, enter the directory ``gnumakefiles``::

    $ cd SALMON/gnumakefiles

In the directory, ``Makefile`` files are prepared for several architectures:

- gnu-mpi
- intel-mpi
- gnu-without-mpi
- intel-without-mpi

``Makefile`` files with ``*-without-mpi`` indicate that they are for single processor environment.
Choose ``Makefile`` appropriate for your environment, and execute the make command::

    $ make -f Makefile.PLATFORM

If the make proceeds successful, a binary file is created in the directory ``SALMON/bin/``.


.. _troubleshooting-install:

Troubleshooting of the Installation Process
-------------------------------------------

.. _installation-cmake:

Installation of CMake
~~~~~~~~~~~~~~~~~~~~~

The `CMake <https://cmake.org/>`_ is a cross-platform build tool. In order to build the
SALMON from the source code, the CMake of version 3.14.0 or later is
required. You may install it following one of the three instructions
below.


Installation by package manager
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your system has a built-in package manager, you may conveniently
install the CMake tools as below:

**Debian/Ubuntu Linux**

::

   sudo apt-get install cmake

**Fedora Linux/CentOS**

::

   sudo yum install cmake

**openSUSE Linux**

::

   sudo zypper install cmake


Installation from source code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can get the source code distribution from the `download page <https://cmake.org/download/>`__. In
this time, we will use the cmake version 3.16.8 as an example. Download
the archive by ``wget`` comamnd and unpack it as below:

::

   wget https://cmake.org/files/v3.16/cmake-3.16.8.tar.gz
   tar -zxvf cmake-3.16.8.tar.gz

And, move to the unpacked directory and build.

::

    
   cd cmake-3.16.8
   ./configure --prefix=INSTALLATION_DIRECTORY
   make
   make install

(replace ``INSTALLATION_DIRECTORY`` to your installation directory.)

Next, to utilize the ``cmake`` command, it is required that the
executable are settled inside the directory specified in your ``$PATH``.
If you use the bash shell, edit ``~/.bashrc`` and append the line:

::

   export PATH=INSTALLATION_DIRECTORY/bin:$PATH

and reload the configuration:

::

   source ~/.bashrc

