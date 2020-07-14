.. _for_large_scale_simulation:

For executing large-scale simulation
====================================

This page will describes the note of large-scale simulation under world top-level supercomputers.
Probably, almost all of users no need to read the page.
If you want the fast simulation on large-scale supercomputers, the page will helps to improve the calculation performance.

MPI process distribution manner
-------------------------------

SALMON provides three variables to determine the process distribution/allocation.

- ``nproc_k``
- ``nproc_ob``
- ``nproc_rgrid(3)``

As the default, SALMON will decides the process distribution automatically.
However, in many case, the manner will not achieves better performance than explicit assignment.

Our recommended chart of the process distribution as below,

If you use the k-point (k-point is greater than 1):

  - First, process assigns to ``nproc_k``.
  - Remaining process assigns to ``nproc_ob``.
  - ``nproc_rgrid = 1, 1, 1`` is best case of this simulation

    - We assume real-space grid (``num_rgrid``) is very small such as 16^3.
 
Else:

  - First, process assigns to ``nproc_ob``.
  - Remaining process assigns to ``nproc_rgrid``.

    - If real-space grid size (``num_rgrid(1:3) = al(1:3) / dl(1:3)``) is roughly larger than 64^3, please you consider the balanced distribution manner between ``nproc_rgrid`` and ``nproc_ob``.

Improve the performance of eigenvalues solver
---------------------------------------------

In DFT calculation, eigenvalues solver is the performance bottleneck in the entire simulation.
We recommend the ScaLAPACK or EigenExa are used to calculate large-scale simulation.
You should rebuild the SALMON with ScaLAPACK/EigenExa enabling, could you please read the :any:`install-and-run`.

When executing SALMON, ``yn_scalapack = 'y'`` or ``yn_eigenexa = 'y'`` should be written to your inputfile::

  &parallel
    yn_scalapack = 'y'               ! use ScaLAPACK to solve eigenvalues
    !yn_eigenexa  = 'y'              ! use EigenExa
    yn_diagonalization_red_mem = 'y' ! to optimize the reducing memory consumption
  /

ScaLAPACK/EigenExa solves eigenvalues problem with ``nproc_ob`` process distribution.
If ``nproc_ob = 1``, ScaLAPACK/EigenExa will performs as same as LAPACK library.

Improve the performance of Hartree solver
-----------------------------------------

We provide the FFT solver to calculate Poisson equation (Hartree potential) in periodic system.
You should be satisfy as following conditions to use FFTE::

  num_rgrid(1) mod nproc_rgrid(2) = 0
  num_rgrid(2) mod nproc_rgrid(2) = 0
  num_rgrid(2) mod nproc_rgrid(3) = 0
  num_rgrid(3) mod nproc_rgrid(3) = 0

  The prime factors for number of real-space grids (num_rgrid(1:3)) must be combination of 2, 3 or 5.


And ``yn_ffte = 'y'`` is set to input file::

  &parallel
    yn_ffte = 'y'
  /

Improve IO performance (write/read wavefunction)
------------------------------------------------

The almost all of supercomputer provides the distributed filesystem such as Lustre.
Distributed filesystem has a meta-data server (MDS) and object-storage server (OST).
OST stores the real user data files, and MDS has the address of OST where user data file.
When accessing to data files in OST, the process send a query of the OST address to MDS.
Therefore, the network contention maybe occurs by the query to MDS.

In many implementation of the filesystem, which MDS to query is determined by the directory structure.
``method_wf_distributor`` and ``nblock_wf_distribute`` are parameters to reduce network contention::

  &control
    method_wf_distributor = 'slice' ! wavefunction divided stores to a file every orbital function.
    nblock_wf_distribute  = 32      ! 32 orbital function files store a same directory.
  /

Improve the communication performance for mesh-torus network system
-------------------------------------------------------------------

Large-scale supercomputer often adopts a mesh-torus network system such as Cray dragon-fly and Fujitsu Tofu due to expensive costs.
We implement the specialized MPI process distribution (communicator creation rule) to improve the performance under large-scale mesh-torus network system.

Currently, we provide the communicator creation rule for "Supercomputer Fugaku", which is developed by RIKEN R-CCS and Fujitsu limited.
Fugaku provides a 6-D mesh-torus network which call as "Tofu-D", user can controls it as 3-D logical network.
SALMON computes 5-D array (wavefunction(x, y, z, orbital, k-point)) domain, therefore, we desire the maps 3-D network to 5-D array distribution.

We defined as following variables and conditons to assign 3-D mesh-torus network to 5-D array distribution::

  PW           = nproc_ob * nproc_k
  (PX, PY, PZ) = nproc_rgrid
  PPN          = '# of process per node' (we recommend value is 4 in Fugaku)
  
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

By these condition, you can determine the suitable process distribution and Tofu-D network shape (compute node shape).
``process_allocation`` input variable controls the order of process distribution.
It improves dominant communication performance when executing multiple processes each node.

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

