.. _install-and-run:

Install and Run
================

Prerequisites
----------------

In this guide, it is assumed that readers have a basic knowledge of Linux and its command line operations.
For the installation of SALMON, following packages are required.

- Fortran90/C compiler. SALMON assumes users have one of the following compilers:

  - GCC (GNU Compiler Collection)
  - Intel Compiler
  - Fujitsu Compiler (at FX100 and A64FX)
  - Nvidia HPC SDK Compiler

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
You can also get the file by::

  $ wget http://salmon-tddft.jp/download/SALMON-<VERSION>.tar.gz

To extract files from the downloaded file ``SALMON-<VERSION>.tar.gz``, type the following command in the command-line::

  $ tar -zxvf ./SALMON-<VERSION>.tar.gz

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

==================  =======================================  ===================  =================
arch                Detail                                   Compiler             Numerical Library
==================  =======================================  ===================  =================
intel-oneapi        Intel oneAPI (cross-architecture)        Intel Compiler       Intel MKL
intel-knl           Intel Knights Landing                    Intel Compiler       Intel MKL
intel-knc           Intel Knights Corner                     Intel Compiler       Intel MKL
intel-avx           Intel Processer (Ivy-, Sandy-Bridge)     Intel Compiler       Intel MKL
intel-avx2          Intel Processer (Haswell, Broadwell ..)  Intel Compiler       Intel MKL
intel-avx512        Intel Processer (Skylake-SP)             Intel Compiler       Intel MKL
fujitsu-fx100       FX100 Supercomputer                      Fujitsu Compiler     SSL-II
fujitsu-a64fx-ea    A64FX processor (Fugaku, FX1000, FX700)  Fujitsu Compiler     SSL-II
nvhpc-openmp        Nvidia OpenMP (CPU)                      Nvidia HPC Compiler  Nvidia HPC SDK
nvhpc-openacc       Nvidia OpenACC (GPU)                     Nvidia HPC Compiler  Nvidia HPC SDK
nvhpc-openacc-cuda  Nvidia OpenACC+CUDA (GPU)                Nvidia HPC Compiler  Nvidia HPC SDK
==================  =======================================  ===================  =================

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

SALMON utilizes norm-conserving (NC) pseudpotentials. 
Filenames of pseudopotentials should be written in the input file.

You may find pseudopotentials of some elements in the samples prepared in :any:`Exercises`.
In SALMON, several formats of pseudopotentials may be usable (listed below).
For example, pseudopotentials with an extension ``.fhi`` can be obtained from the ABINIT website (this is a part of previous atomic data files for the ABINIT code).

=========================================================  =============  =====================================================================================
Pseudopotential                                            extension      Website
=========================================================  =============  =====================================================================================
Fritz-Haber-Institute (FHI) pseudopotentials               ``.fhi``       https://abinit.github.io/abinit_web/ATOMICDATA/LDA_FHI.zip 
                                                                          (for LDA), 
                                                                          https://abinit.github.io/abinit_web/ATOMICDATA/fhi.zip
                                                                          (for GGA) 
Pseudopotentials for the OpenMX code                       ``.vps``       https://www.openmx-square.org/vps_pao2019/
Format 8 for ABINIT norm-conserving pseudopotentials       ``.psp8``      https://abinit.github.io/abinit_web/pseudopotential.html , 
                                                                          http://www.pseudo-dojo.org/
Unified-pseudopotential-format (NC type only in SALMON)    ``.upf``       http://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format , 
                                                                          http://www.pseudo-dojo.org/
=========================================================  =============  =====================================================================================


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
- Executable file is prepared as ``salmon`` in the standard build procedure.
- To start calculations, ``inputfile.inp`` should be read through ``stdin``.


MPI process distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SALMON provides three variables to determine the process distribution/allocation.

- ``nproc_k``
- ``nproc_ob``
- ``nproc_rgrid(3)``

In SALMON, the process distribution is determined automatically as default.
However, in many situations, an explicit assignment of the process distribution
will provide a better performance than the default setting.

We recommend to distribute the processes as follows,

If you use k-points ( the number of k-points is greater than 1) and the number of 
the real-space grid (``num_rgrid``) is not very large (about 16^3):

  - First, assign many processes to ``nproc_k``.
  - Then, assign the remaining processes to ``nproc_ob``.
  - Not dividing the spatial grid,  ``nproc_rgrid = 1, 1, 1``.
 
Else:

  - First, assign the processes to ``nproc_ob``.
  - Then, assign the remaining processes to ``nproc_rgrid``.

    - If real-space grid size (``num_rgrid(1:3) = al(1:3) / dl(1:3)``) is equal to or larger than about 64^3, 
    you should find a balanced distribution between ``nproc_rgrid`` and ``nproc_ob``.


.. _for_large_scale_simulation:

Tips for large-scale calculation
-----------------------------------

We explain below some tips that will be useful to improve performance when you carry out 
large scale simulations using world top-level supercomputers.
Therefore, the following contents will only be useful only for limited users.

Improve the performance of the eigenvalues solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In DFT calculations of large systems, subspace diagonalization becomes the performance bottleneck
in the entire calculation. Therefore, it is important to use a parallel eigenvalues solver.
In SALMON, a LAPACK routine without parallelization is used for the diagonalization as default.
As parallelized solvers, ScaLAPACK and EigenExa are usable.
To use them, it is necessary to rebuild SALMON enabling ScaLAPACK/EigenExa.
You can find the instruction in :any:`additional-options-in-configure`.

To execute SALMON using ScaLAPACK/EigenExa, either ``yn_scalapack = 'y'`` or ``yn_eigenexa = 'y'`` should be 
included in the inputfile::

  &parallel
    yn_scalapack = 'y'         ! use ScaLAPACK for diagonalization
    !yn_eigenexa  = 'y'        ! use EigenExa
  /

ScaLAPACK/EigenExa solves the eigenvalue problem with ``nproc_ob`` process distribution.
If ``nproc_ob = 1``, ScaLAPACK/EigenExa will perform in the same way as the LAPACK library.

Improve the performance of Hartree solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For periodic systems, a Fourier transformation is used to solve the Poisson equation (to calculate the Hartree potential).
In SALMON, a simple Fourier transformation without Fast Fourier Transformation (FFT) is used as default.
In SALMON, a parallelized FFT routine, FFTE, is usable and works efficiently for large systems.
In using FFTE, the following conditions should be satisfied::

  num_rgrid(1) mod nproc_rgrid(2) = 0
  num_rgrid(2) mod nproc_rgrid(2) = 0
  num_rgrid(2) mod nproc_rgrid(3) = 0
  num_rgrid(3) mod nproc_rgrid(3) = 0

  In addition, the prime factors for the number of real-space grid of each direction (num_rgrid(1:3)) must be a combination of 2, 3 or 5.


To use FFTE, ``yn_ffte = 'y'`` should be included in the input file::

  &parallel
    yn_ffte = 'y'
  /

Improve IO performance (write/read wavefunction)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Almost all supercomputer systems provide distributed filesystems such as Lustre.
Distributed filesystems are equipped with a meta-data server (MDS) and an object-storage server (OST).
The OST stores real user data files, and the MDS stores the address of the user date files in the OST.
When accessing to the data files in the OST, the process send a query about the OST address to MDS.
Then, a network contention may occur in the query process.

In most implementations of the filesystem, the MDS that replies to the query is determined by the directory structure.
For a calculation in which k-point is not used, 
``method_wf_distributor`` and ``nblock_wf_distribute`` are prepared to reduce the network contention::

  &control
    method_wf_distributor = 'slice' ! every orbital function is stored as a single file.
    nblock_wf_distribute  = 32      ! files of 32 orbital functions are stored in one directory.
  /

Improve the communication performance for mesh-torus network system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Large-scale supercomputers often adopt a mesh-torus network system such as Cray dragon-fly and Fujitsu Tofu to achieve
high scalability with relatively low cost. 
In SALMON, a special MPI process distribution (communicator creation rule) is prepared to improve the performance 
in large-scale mesh-torus network systems.

Currently, we provide the communicator creation rule for "Supercomputer Fugaku", 
which is developed by RIKEN R-CCS and Fujitsu limited.
Fugaku is equipped with a 6-D mesh-torus network which is called "Tofu-D". 
Users may control it as a 3-D logical network.
SALMON utilizes 5-D array (wavefunction(x, y, z, orbital, k-point)) as a domain for parallelization.
We create a map that connects the 3-D network to the 5-D array distribution.

We introduce the following variables and conditons to assign the 3-D mesh-torus network to the 5-D array distribution::

  PW           = nproc_ob * nproc_k
  (PX, PY, PZ) = nproc_rgrid
  PPN          = '# of process per node' (we recommend the value 4 in Fugaku)
  
  Requested process shape: (PX, PY, PZ, PW)
  Tofu-D network    shape: (TX, TY, TZ)
  Actual process    shape: (TX * PPN, TY, TZ)

  if (process_allocation == 'grid_sequential'):
    PW  = PW1 * PW2 * PW3
    PW1 = (TX * PPN) / PX
    PW2 = TY         / PY
    PW3 = TZ         / PZ
    TX  = (PX * PW1) / PPN
    TY  = PY * PW2
    TZ  = PZ * PW3

  else if (process_allocation == 'orbital_sequential'):
    PX  = PX1 * PX2 * PX3
    PX1 = (TX * PPN) / PW
    PX2 = TY         / PY
    PX3 = TZ         / PZ
    TX  = (PW * PX1) / PPN
    TY  = PY * PX2
    TZ  = PZ * PX3

From these conditions, you can determine the suitable process distribution and the Tofu-D network shape (compute node shape).
``process_allocation`` input variable controls the order of the process distribution.
It indicates which communications should be executed in closer processes.

- ``process_allocation = 'grid_sequential'``

  - ``(PX, PY, PZ, PW)``, ``nproc_rgrid`` major ordering
  - improves ``nproc_rgrid`` related communication performance
  - communicator: ``s_parallel_info::icomm_r, icomm_x, icomm_y, icomm_z, icomm_xy``
  - suitable ``theory``: ``'dft'`` and ``'dft_md'``

- ``process_allocation = 'orbital_sequential'``

  - ``(PW, PY, PZ, PX)``, ``nproc_ob`` major ordering
  - improves ``nproc_ob`` related communication performance
  - communicator: ``s_parallel_info::icomm_o and icomm_ko``
  - suitable ``theory``: ``'tddft_response', 'tddft_pulse', 'single_scale_maxwell_tddft'`` and ``'multi_scale_maxwell_tddft'``

.. _GPU:

GPU acceleration
~~~~~~~~~~~~~~~~~~~~~

GPU acceleration (OpenACC or OpenACC+CUDA) for the DFT/TDDFT computation is available. 
For compiling SALMON for GPUs, specify ``--arch=nvhpc-openacc`` (OpenACC, recommended) or ``--arch=nvhpc-openacc-cuda`` (OpenACC+CUDA) option when executing ``configure.py``.
This option is currently under development and tested only for NVIDIA HPC SDK compiler with NVIDIA GPUs.

Multi-GPU run
^^^^^^^^^^^^^^^^^^^^^^^^^

For MPI calculations with multiple GPUs, the assignment of MPI processes to GPUs via CUDA_VISIBLE_DEVICES and the use of nvidia-cuda-mps-control can improve the performance of SALMON. The following example is a wrapper script for that::

    $ cat wrapper.sh
    #! /bin/bash
    ### wrapper.sh
    NCUDA_GPUS=${NCUDA_GPUS:-`nvidia-smi -L | wc -l`}
    if $OMPI_COMM_WORLD_LOCAL_SIZE -gt $NCUDA_GPUS 
    then
      if $OMPI_COMM_WORLD_LOCAL_RANK -eq 0  
      then
        nvidia-cuda-mps-control -d
      fi
      sleep 10
    fi
    export CUDA_VISIBLE_DEVICES=$((${OMPI_COMM_WORLD_LOCAL_RANK} % ${NCUDA_GPUS})) 
    exec $@
    if $OMPI_COMM_WORLD_LOCAL_SIZE -gt $NCUDA_GPUS 
    then
      echo quit | nvidia-cuda-mps-control 
    fi

Here, we used environment variables of OpenMPI, such as $OMPI_COMM_WORLD_LOCAL_SIZE.
For MPI execution, use the following command::

    $ mpirun -np ${num_MPI_processes} -npernode ${num_MPI_processes_per_node} \
       wrapper.sh ${program} < ${input} |& tee log.stdout

Here, ${program} is the path of SALMON, ${input} is the input file, etc.

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


Appendix
------------

.. _additional-options-in-configure:

Additional options in configure.py script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Manual specifications of compiler and environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When executing ``configure.py``, you can manually specify the compiler and environment variables instead of specifying the architecture ``--arch``, for example::

    $ python ../configure.py FC=mpiifort CC=mpiicc FFLAGS="-xAVX" CFLAGS="-restrict -xAVX" --enable-mpi

The list of options of ``configure.py`` can be found by::

    $ python ../configure.py --help

The major options are as follows:

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
LDFLAGS                                  linker flags
=======================================  ===================================================

Required libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the build procedure of SALMON, CMake searches the following libraries.
If the libraries are not found in the path specified by environment variables, the required libraries will be downloaded and compiled automatically.

- BLAS/LAPACK 

    - Required by default compilation.
    - Most math libraries include BLAS/LAPACK by default.
    - ``--with-lapack``: Path specification.
    - If the library is not found, it will be automatically downloaded from http://www.netlib.org/lapack/

- ScaLAPACK

    - Required by ``--enable-scalapack``.
    - ``--with-lapack``: Path specification.
    - If the library is not found, it will be automatically downloaded from http://www.netlib.org/scalapack/

- Libxc

    - Required by ``--enable-libxc``.
    - ``--with-libxc``: Path specification.
    - If the path is unspecified, the library will be automatically downloaded from https://libxc.gitlab.io

- EigenExa

    - Required by ``--enable-eigenexa``. (``--enable-scalapack`` is also required for EigenExa.)
    - EigenExa will be downloaded and built automatically even if the library is installed on your machine.
    - Automatically download from https://www.r-ccs.riken.jp/labs/lpnctrt/assets/img/EigenExa-2.4b.tgz

Build for single process calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the ``--arch`` option, MPI parallelization is enabled as default.
If you use a single processor machine, explicitly specify ``--disable-mpi`` in executing the python script::

    $ python ../configure.py --arch=<ARCHITECTURE> --disable-mpi


Build by GNU Compiler Collection (GCC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The architecture option ``--arch`` does not support GNU Compiler Collection (GCC).
If you want to build SALMON by GCC, specify ``FC`` and ``CC`` flags as follows::

    $ python ../configure.py FC=gfortran CC=gcc --enable-mpi

Here, ``--enable-mpi`` is required for the MPI parallelization.
Note that the MPI parallelization is disabled as default when ``--arch`` option is not used.
Compiler options can also be specified by ``FFLAGS`` and ``CFLAGS``. For GCC 10 or later versions, ``FFLAGS="-fallow-argument-mismatch"`` may be required.


Compilation examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some compilation (configure) examples in several environments are shown below.

- Wisteria-Odyssey (University of Tokyo) & Fujitsu compiler, compiling with EigenExa::

    $ python3 ../configure.py --arch=fujitsu-a64fx-ea --enable-scalapack --enable-eigenexa FFLAGS="-fPIC"

- Cygnus (GPU supercomputer @ University of Tsukuba) & NVidia HPC SDK compiler version 23.11::

    $ module load cmake/3.18.6 openmpi/nvhpc/23.11
    $ python3 ../configure.py --arch=nvhpc-openacc LDFLAGS=-L/system/apps/nvhpc/23.11/Linux_x86_64/23.11/math_libs/lib64/

- AWS Graviton2 machine (Amazon EC2 T4g instance) & Arm compiler::

    $ python3 ../configure.py FC=armflang CC=armclang FFLAGS="-armpl" CFLAGS="-armpl"

- MacOS & GCC version 11::

    $ brew install gcc@11
    $ export FC=/opt/homebrew/Cellar/gcc@11/11.5.0/bin/gfortran-11
    $ export CC=/opt/homebrew/Cellar/gcc@11/11.5.0/bin/gcc-11
    $ export CXX=/opt/homebrew/Cellar/gcc@11/11.5.0/bin/g++-11
    $ python ../configure.py FFLAGS="-fallow-argument-mismatch"

.. _FFTW:

Compilation with FFTW library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For solving the Poisson equation for the Hartree potential, SALMON uses the discrete Fourier transform. 
FFTW library (https://www.fftw.org) is available for fast calculation. 
When executing ``configure.py``, specify ``--enable-fftw`` option and linker flags for FFTW such as ``LDFLAGS="-lfftw3_mpi -lfftw3"``.

Exapmle::

    $ python ../configure.py --arch=ARCHITECTURE --enable-fftw LDFLAGS="-lfftw3_mpi -lfftw3"

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


