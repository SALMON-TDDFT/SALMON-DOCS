.. _List of input keywords:

List of input keywords
======================


'[Trial]' : These options are not tested well

.. _&calculation:

&calculation
------------

.. _theory:

theory
^^^^^^


character, default=''

   | Choice of a theory to be used in the calculation.
   | Options:
   |   ``dft``  / ground state calculation based on DFT
   |   ``dft_md``  / ab initio MD simulations based on DFT (electronic ground state)
   |   ``tddft_response``  / linear response TDDFT calculation in real time
   |   ``tddft_pulse``  / simulations under pulsed electric field based on TDDFT
   |   ``single_scale_maxwell_tddft``  / single-scale simulation coupling Maxwell and TDDFT
   |   ``multi_scale_maxwell_tddft``  / multiscale simulation coupling Maxwell and TDDFT
   |   ``maxwell``  / electromagnetic analysis using finite difference time domain (FDTD) method
   |   ``dft_k_expand`` / convert checkpoint data of dft with k-points calculation to that of larger supercell system with gamma-point
   |   ``sbe`` / [Trial] simulations under pulsed electric field based in semiconductor Bloch equation 
   |   ``maxwell_sbe`` / [Trial] multiscale simulation coupling Maxwell and semiconductor Bloch equation

.. _yn_md:

yn_md
^^^^^

[Trial] character, default='n'

   | Available for ``theory='dft'`` (ground-state MD) and ``theory='tddft_pulse'`` (Ehrenfest MD).
   | Switch for molecular dynamics calculation.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_opt:

yn_opt
^^^^^^


[Trial] character, default='n'

   | Available for ``theory='dft'``.
   | Switch for geometry optimization.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _&control:

&control
--------

.. _sysname:

sysname
^^^^^^^

character, default='default'

   | Available for all options of ``theory``.
   | A prefix of output files.

.. _base_directory:

base_directory
^^^^^^^^^^^^^^

character, default='./'

   | Available for all options of ``theory``.
   | Name of a directory where major output files are stored.

.. _yn_restart:

yn_restart
^^^^^^^^^^

character, default='n'

   | Available for the DFT/TDDFT based options of ``theory`` and ``theory='maxwell'``.
   | Whether to continue previous calculation (restart) or start a new calculation.
   | Options:
   |   ``'y'`` / enable (restart)
   |   ``'n'`` / disable (new calculation)

.. _directory_read_data:

directory_read_data
^^^^^^^^^^^^^^^^^^^

character, default='restart/'

   | Available for ``yn_restart='y'``.
   | Directory name to read data that are required in the present calculation (restart) and were generated in previous calculations. For TDDFT based options, it specifies the name of the directory containing ground state results that were stored in 'data_for_restart'. When restarting from a checkpoint, it specifies the name of the directory that contains the checkpoint data.

.. _yn_self_checkpoinnt:

yn_self_checkpoint
^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for the DFT/TDDFT based options of ``theory``.
   | With this option, each process writes/reads the restart (and checkpoint) data independently (self data format) so that the restart cost is reduced for large systems. Note that the number of processes and their assignments must be unchanged in restarting. The data is written out into 'checkpoint_gs_XXXXXX/' (DFT) or 'checkpoint_rt_XXXXXX/' (TDDFT).
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _checkpoint_interval:

checkpoint_interval
^^^^^^^^^^^^^^^^^^^

integer, default=-1

   | Available for the DFT/TDDFT based options of ``theory`` and ``theory='maxwell'``.
   | Interval of time steps (iteration steps) to write down the checkpoint data during the time-propagation (SCF iteration). Checkpoint data will not be written if a negative value is set.

.. _yn_reset_step_restart:

yn_reset_step_restart
^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``yn_restart='y'`` in the DFT/TDDFT based options of ``theory``.
   | With this option, the counter of the SCF iteration step (for DFT) or the counter of the time propagation step (for TDDFT) is reset to 0 at the restart. In the SCF iteration, the density data in the previous SCF iteration step are abondoned.

.. _read_gs_restart_data:

read_gs_restart_data
^^^^^^^^^^^^^^^^^^^^

character, default='all'

   | Available for ``yn_restart='y'`` with ``theory='dft'``.
   | Specify which data are read in the restart. Specified data that are generated in the previous calculation and are contained in the restart (or checkpoint) directory are used in restarting the SCF iteration of DFT. The default option ``'all'`` indicates the complete restart. In other options, a part of restart data are used (other data are prepared in the same way as in the initial SCF step).
   | Options:
   |   ``all``  / all of restart data are read
   |   ``all:single``  / same as ``all`` option but the data is read in the single file format even though the self data format is specified with ``yn_self_checkpoint='y'`` (i.e., the restart data is read in the single file format while written out in the self format)
   |   ``rho_inout``  / only electron densities including those of previous iteration steps are read (from rho_inout.bin file)
   |   ``rho_inout:single``  / same as ``rho_inout`` option but the data is read in the single file format even though the self data format is specified with ``yn_self_checkpoint='y'``
   |   ``rho``  / only the latest electron density is read (from user-made data)
   |   ``wfn``  / only orbital wavefunctions are read

.. _write_gs_restart_data:

write_gs_restart_data
^^^^^^^^^^^^^^^^^^^^^

character, default='all'

   | Available for ``theory='dft'``.
   | Options
   |   ``all``  / all of restart data are written out
   |   ``rho_inout``  / only electron densities including those of previous iteration steps are written out
   |   ``wfn``  / only orbital wavefunctions are written out
   |   ``checkpoint_only`` / the restart data are outputted only in the self data format (separated data for each process) at the last step into 'checkpoint_gs_XXXXXX/' directory (``yn_self_checkpoint='y'`` is required) without generating the restart data into 'data_for_restart/' directory in the single file format.
   | Output data files are written out in the restart (or checkpoint) directory.
   | The default option ``'all'`` gives the complete set of restart data.

.. _time_shutdown:

time_shutodown
^^^^^^^^^^^^^^

[Trial] real(8), default=-1d0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Timer for automatic shutdown. The unit is second.
   | If a negative time is set, the automatic shutdown will not be performed.

.. _method_wf_distributor:

method_wf_distributor
^^^^^^^^^^^^^^^^^^^^^

character, default='single'

   | Available for the DFT/TDDFT based options of ``theory``.
   | A method to save/load orbital wavefunctions.
   | Options
   |   ``single``: all orbital wavefunctions are saved(loaded) to(from) a single file.
   |   ``slice`` : each orbital wavefunction is saved(loaded) to(from) a file. This choice reduces the I/O costs, and increase the flexiblility in handling files for large systems.

.. _nblock_wf_distribute:

nblock_wf_distribute
^^^^^^^^^^^^^^^^^^^^

integer, default='16'

   | Available for ``method_wf_distributor='slice'``.
   | In the 'slice' mode, ``nblock_wf_distribute`` files are saved in one directory.

.. _&units:

&units
------

.. _unit_system:

unit_system
^^^^^^^^^^^

character, default='au'

   | Unit system to be used in input variables and some of output files.
   | If ``unit_system = 'A_eV_fs'`` is chosen, Angstrom for length, eV for energy, and fs for time are adopted.
   | For isolated systems specified by ``yn_periodic = 'n'`` in ``&system``, a unit of 1/eV is used for the output files of DOS and PDOS if ``unit_system = 'A_eV_fs'`` is specified, while atomic unit is used if not. For other output files, the Angstrom/eV/fs units are used irrespective of the input keyword. For periodic systems specified by ``yn_periodic = 'n'`` in ``&system``, the unit system specified by this input keyword is used for most output files. To confirm the unit, see the first few lines of output files.
   | Options:
   |   ``'au'`` or ``'a.u.'`` / atomic unit system
   |   ``'A_eV_fs'`` / Angstrom-eV-fs unit system

.. _&parallel:

&parallel
---------

.. _nproc_k:

nproc_k
^^^^^^^

.. _nproc_ob:

nproc_ob
^^^^^^^^

.. _nproc_rgrid(3):

nproc_rgrid(3)
^^^^^^^^^^^^^^

integer, default=0

   | Options:
   |   ``nproc_k``/ Number of MPI parallelization for k-points of electron orbitals.
   |   ``nproc_ob``/ Number of MPI parallelization for orbital index of electron orbitals.
   |   ``nproc_rgrid(3)'``/ Number of MPI parallelization for each direction of real-space grid that are used for electron orbitals and density.
   |
   | Defaults are ``0`` for ``nproc_k``/``nproc_ob`` and ``(0,0,0)`` for ``nproc_rgrid``. In the default choice, MPI assignment is achieved atomatically. Users can specify ``nproc_k``, ``nproc_ob``, and ``nproc_rgrid`` manually. In that case, there are several constraints that should be fulfilled:
   |   ``nproc_k`` must be set to ``1`` for ``&system/yn_periodic='n'``.
   |   ``nproc_k`` and ``nproc_ob`` must be set to ``1`` for ``theory='maxwell'``.
   |   ``nproc_k`` \* ``nproc_ob`` \* ``nproc_rgrid(1)`` \* ``nproc_rgrid(2)`` \* ``nproc_rgrid(3)`` \= total number of processes.

.. _yn_ffte:

yn_ffte
^^^^^^^

character, default='n'

   | Available for the DFT/TDDFT based options of ``theory`` with ``&system/yn_periodic='y'``
   | For periodic systems, SALMON uses Fourier transformation to solve a poisson equation.
   | This switch selects if FFTE library is used or not. If FFTE is not used, the Fourier transformation in a simple algorithm is carried out.
   | Options
   |   ``'y'`` / enable
   |   ``'n'`` / disable
   |
   | To enable it, following relations must be satisfied.
   |   ``mod(num_rgrid(1), nproc_rgrid(2)) == 0``
   |   ``mod(num_rgrid(2), nproc_rgrid(2)) == 0``
   |   ``mod(num_rgrid(2), nproc_rgrid(3)) == 0``
   |   ``mod(num_rgrid(3), nproc_rgrid(3)) == 0``

.. _yn_fftw:

yn_fftw
^^^^^^^

character, default='n'

   | Available for the DFT/TDDFT based options of ``theory`` with both ``&system/yn_periodic='y'`` and ``&system/yn_periodic='n'``.
   | For isolated systems, this option is effective when ``&poisson/method_poisson='ft'``
   | This switch selects if FFTW library is used or not. If FFTW is not used, the discrete Fourier transformation in a simple algorithm is carried out.
   | Caution: This variable is effective only when ``--enable-fftw`` is specified at the configure.
   | Options
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_scalapack:

yn_scalapack
^^^^^^^^^^^^

character, default='n'

   | Available for ``&calculation/theory='dft' or 'dft_md'``
   | To calculate large systems, an eigenvalue problem in the subspace diagonalization becomes a bottle-neck in the ground state calculation. In SALMON, ScaLAPACK library can be used to solve the eigenvalue problem.
   | To enable it, it is necessary to link ScaLAPACK library when you build SALMON.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_gramschmidt_blas:

yn_gramschmidt_blas
^^^^^^^^^^^^^^^^^^^

character, default='y'

   | Available for ``&calculation/theory='dft' or 'dft_md'``
   | This switch selects if BLAS library is used or not in Gram Schmidt routines.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_eigenexa:

yn_eigenexa
^^^^^^^^^^^

character, default='n'

   | Available for ``&calculation/theory='dft' or 'dft_md'``
   | SALMON can use RIKEN R-CCS EigenExa library to solve eigenvalue problem in subspace diagonalization. It is more efficient than ScaLAPACK to diagonalize matrices of large dimension. To enable it, it is necessary to link both ScaLAPACK and EigenExa libraries when you build SALMON.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_diagonalization_red_mem:

yn_diagonalization_red_mem
^^^^^^^^^^^^^^^^^^^^^^^^^^

character, Default='n'

   | Available for ``&parallel/yn_scalapack='y'`` or ``&parallel/yn_eigenexa='y'``
   | This option reduces memory consumption in using ScaLAPACK/EigenExa libraries.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _process_allocation:

process_allocation
^^^^^^^^^^^^^^^^^^

character, default='grid_sequential'

   | This controlls the order of process allocation.
   | Options:
   |   ``'grid_sequential'``    / real-space grid major ordering.
   |   ``'orbital_sequential'`` / orbital-space major ordering.
   |
   | Suggestion:
   |   ``&calculation/theory='dft' or 'dft_md'``            / ``'orbital_sequential'``
   |   ``&calculation/theory='tddft*' or '*maxwell_tddft'`` / ``'grid_sequential'``

.. _&system:

&system
-------

.. _yn_periodic:

yn_periodic
^^^^^^^^^^^

character, default='n'

   | Available for all options of ``theory``.
   | Specify boundary condition for electron orbitals.
   | Options:
   |   ``'y'`` / periodic systems (crystalline solids)
   |   ``'n'`` / isolated systems (molecules and nano-particles)

.. _spin:

spin
^^^^

character, default='unpolarized'

   | Available for the DFT/TDDFT based options of ``theory``.
   | It specifies the spin state of the system, spin-unpolarized (closed shell) or spin-polarized (open shell).
   | Options
   |   ``'unpolarized'`` / spin-unpolarized systems (default)
   |   ``'polarized'`` / spin-polarized systems
   |   ``'noncollinear'`` / noncollinear spin systems (see ``yn_spinorbit``)

.. _al(3):

al(3)
^^^^^

real(8), default=0d0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Spatial box size or lattice constants for cuboid cell (x, y, z).
   | For nonorthogonal cell, see ``al_vec1(3)``, ``al_vec2(3)``, ``al_vec3(3)``.

.. _al_vec1(3):

al_vec1(3)
^^^^^^^^^^

.. _al_vec2(3):

al_vec2(3)
^^^^^^^^^^

.. _al_vec3(3):

al_vec3(3)
^^^^^^^^^^

real(8), default=0d0

   | Available for ``yn_periodic = 'y'`` in the DFT/TDDFT based options of ``theory``.
   | Primitive lattice vectors for nonorthogonal cell. For cuboid cell, see ``al(3)``.

.. _nstate:

nstate
^^^^^^

integer, default=0

   | Available for the DFT/TDDFT based options of ``theory``.
   |  of orbitals/bands to be calculated. In the time evolution calculation of dielectrics, only occupied orbitals are evolved even when more ``nstate`` is specified.

.. _nelec:

nelec
^^^^^

integer, default=0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Number of valence electrons in the system.

.. _nelec_spin(2):

nelec_spin(2)
^^^^^^^^^^^^^

integer, Default=0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Number of up/down-spin electrons are specified by ``nelec_spin(1)/nelec_spin(2)``.
   | This option is incompatible with ``nelec``. (If ``nelec_spin`` is specified, ``nelec`` is ignored.)

.. _temperature:

temperature
^^^^^^^^^^^

real(8), default=-1d0

   | Available for DFT-based options of ``theory``.
   | It specifies the temperature for electrons. The value must be given using the unit of energy as specified in ``&units/unit_system``.
   | The kelvin unit can also be used by the keyword ``temperature_k`` instead of ``temperature`` (see next).
   | Occupation numbers are updated in every SCF step in the following way.
   |   ``temperature < 0`` / the occupation numbers are fixed by ``nelec`` (appropriate for systems with energy gap).
   |   ``temperature = 0`` / redistribution of the occupation numbers by the step function (metallic system at zero temperature).
   |   ``temperature > 0`` / redistribution of the occupation numbers by the Fermi-Dirac distribution function.

.. _temperature_k:

temperature_k
^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for DFT-based options of ``theory``.
   | The same as ``temperature`` but kelvin is used as the unit.

.. _nelem:

nelem
^^^^^

integer, default=0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Number of atomic elements in the system.

.. _natom:

natom
^^^^^

integer, default=0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Number of atoms in the system.

.. _file_atom_red_coor:

file_atom_red_coor
^^^^^^^^^^^^^^^^^^

[Trial] character, default='none'

   | Available for the DFT/TDDFT based options of ``theory``.
   | Name of the file that contains atomic positions given in reduced coordinates. This option is incompatible with ``&system/file_atom_coor``, ``&atomic_coor``, and ``&atomic_red_coor``.

.. _file_atom_coor:

file_atom_coor
^^^^^^^^^^^^^^

[Trial] character, default='none'

   | Available for the DFT/TDDFT based options of ``theory``.
   | Name of the file that contains atomic Cartesian coordinates (The unit is specified by ``&units/unit_system``). This option is incompatible with ``&system/file_atom_coor``, ``&atomic_coor``, and ``&atomic_red_coor``.

.. _yn_spinorbit:

yn_spinorbit
^^^^^^^^^^^^

character, default='n'

   | Available for the DFT/TDDFT based options of ``theory``.
   | Option for the spin-orbit coupling using the j-dependent pseudopotential formalism [Theurich & Hill, PRB 64, 073106 (2001)]. For pseudopotential(s), the UPF or VPS file format is required.

   | Options
   |   ``'y'`` / enable (``spin='noncollinear'`` is required. For ``theory='dft’`` mode, ``method_mixing='simple’`` is recommended.)
   |   ``'n'`` / disable (default)

.. _yn_symmetry:

yn_symmetry
^^^^^^^^^^^

[Trial] character, default='n'

   | Available for orthogonal cell system with the DFT/TDDFT based options of ``theory``.
   | Symmetry option. Pre-generated input file, "sym.dat", is necessary. (details are not explained in the current manual)

   | Options
   |   (e.g.) ``'yyn'`` / symmetry option is applied for the x and y direction (under applied electric field in the z-direction)
   |   ``'n'`` / disable

.. _absorbing_boundary:

absorbing_boundary
^^^^^^^^^^^^^^^^^^

[Trial] character, default='none'

   | Available for the TDDFT based option of ``theory`` with orthogonal unit cell.
   | Absorbing boundary condition for electrons. (T. Nakatsukasa et al., J. Chem. Phys., 114, 2550 (2001))
   | Options:
   |   ``'none'`` / disable (default)
   |   ``'z'`` / absorbing boundary region is set in z direction for ``'yn_periodic = 'y'``

.. _imagonary_potential_w0:

imaginary_potential_w0
^^^^^^^^^^^^^^^^^^^^^^

real(8), default='0d0'

   | Available when ``absorbing_boundary`` options is not ``'none'``.
   | Strength of the absorbing (imaginary) potential.

.. _imaginary_potential_dr:

imaginary_potential_dr
^^^^^^^^^^^^^^^^^^^^^^

real(8), default='0d0'

   | Available when ``absorbing_boundary`` options is not ``'none'``.
   | Thickness of the absorbing (imaginary) potential. For ``absorbing_boundary='z'``, the absorbing region is 0 < z < ``imagnary_potential_dr`` and ``al(3)``-``imagnary_potential_dr`` < z < ``al(3)``

.. _&atomic_red_coor:

&atomic_red_coor
----------------

   | Atomic coordinates in periodic systems (``'yn_periodoc = 'y'``) are specified in reduced coordinates using the following format:
   |
   |    'Si'	 0.00  0.00  0.00  1
   |    'Si'	 0.25  0.25  0.25  1
   |    ...
   |
   | Here, the information of atoms is ordered in row, the first row for the first atom, etc. The number of rows must be equal to ``&system/natom``. Atomic spicies are written in the first column although they are not used in the calculation. The second, third and fourth columns are reduced coordinates for the first, second and third directions, respectively. The fifth column is a serial number of the atom spieces, which is defined in ``&pseudo``.
   | This option is incompatible with ``&system/file_atom_red_coor``, ``&system/file_atom_coor``, and ``&atomic_coor``.

.. _&atomic_coor:

&atomic_coor
------------

   | Atomic coordinates are specified in the same way as ``atomic_red_coor`` but with length dimension. The unit chosen by ``&units/unit_length`` is applied.
   | This option is incompatible with ``&system/file_atom_red_coor``, ``&system/file_atom_coor``, and ``&atomic_red_coor``.

.. _&pseudo:

&pseudo
-------

.. _izatom(:):

izatom(:)
^^^^^^^^^

integer, default=-1

   | Available for the DFT/TDDFT based options of ``theory``.
   | Atomic number of the element. The size of array is equal to ``&system/nelem``.

.. _file_pseudo(:):

file_pseudo(:)
^^^^^^^^^^^^^^

character, default='none'

   | Available for the DFT/TDDFT based options of ``theory``.
   | File name of the pseudopotential file. The size of array is equal to ``&system/nelem``.

.. _lmax_ps(:):

lmax_ps(:)
^^^^^^^^^^

integer, default=-1

   | Available for the DFT/TDDFT based options of ``theory``.
   | Maximum angular momentum of pseudopotential projectors.
   | If not given, values specified in the pseudopotential file will be used. The size of array is equal to ``&system/nelem``.

.. _lloc_ps(:):

lloc_ps(:)
^^^^^^^^^^

integer, default=-1

   | Available for the DFT/TDDFT based options of ``theory``.
   | Angular momentum of the pseudopotential that will be treated as local. The size of array is equal to ``&system/nelem``.

.. _yn_psmask(:):

yn_psmask(:)
^^^^^^^^^^^^

[Trial] character, default='n'

   | Available for the DFT/TDDFT based options of ``theory``.
   | Fourier filtering for pseudopotentials. The size of array is equal to ``&system/nelem``.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _alpha_mask(:):

alpha_mask(:)
^^^^^^^^^^^^^

[Trial] real(8), default=0.8d0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Parameter for the Fourier filtering of the pseudopotential. The size of array is equal to ``&system/nelem``.

.. _gamma_mask(:):

gamma_mask(:)
^^^^^^^^^^^^^

[Trial] real(8), default=1.8d0)

   | Available for the DFT/TDDFT based options of ``theory``.
   | Parameter for the Fourier filtering of the pseudopotential. The size of array is equal to ``&system/nelem``.

.. _eta_mask(:):

eta_mask(:)
^^^^^^^^^^^

[Trial] real(8), default=15.0d0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Parameter for the Fourier filtering of the pseudopotential. The size of array is equal to ``&system/nelem``.

.. _&functional:

&functional
-----------

.. _xc:

xc
^^

character, default='none'

   | Available for the DFT/TDDFT based options of ``theory``.
   | Exchange-correlation functional to be used.
   | In the present version, functionals 'PZ', 'PZM' and 'TBmBJ' are available for both ``yn_periodic = 'y' and 'n'`` calculations in the adiabatic approximation.
   | Options:
   |   ``'PZ'``: Perdew-Zunger LDA :Phys. Rev. B 23, 5048 (1981).
   |   ``'PZM'``: Perdew-Zunger LDA with modification to improve sooth connection between high density form and low density one. :J. P. Perdew and Alex Zunger, Phys. Rev. B 23, 5048 (1981).
   |   ``'TBmBJ'``: Tran-Blaha meta-GGA exchange with Perdew-Wang correlation. :Fabien Tran and Peter Blaha, Phys. Rev. Lett. 102, 226401 (2008). John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992). This potential is known to provide a reasonable description for the bandgap of various insulators. For this choice, the additional mixing parameter 'cval' may be specified. See below.

.. _cval:

cval
^^^^

real(8), default=-1d0

   | Available for ``xc='TBmBJ'``.
   | Mixing parameter in Tran-Blaha meta-GGA exchange potential. If ``cval`` is set to a minus value, the mixing-parameter is evaluated by the formula in the original paper [Phys. Rev. Lett. 102, 226401 (2008)], :math:`\left\langle |\nabla \rho(\mathbf{r};t)| / \rho(\mathbf{r};t) \right\rangle`. However, note that the value may be different from that in all electron calculations.

.. _cname:

cname
^^^^^

.. _xname:

xname
^^^^^

character, default='none'

   | Available for ``theory='XXX'``.
   | XXX

.. _alibxc:

alibxc
^^^^^^

.. _alibx:

alibx
^^^^^

.. _alibc:

alibc
^^^^^

character, default='none'

   | Available for the DFT/TDDFT based options of ``theory``.
   | Since version 1.1.0, exchange-correlation functionals in Libxc library (http://www.tddft.org/programs/libxc/) have been usable in SALMON. At present, usable functionals are limited to LDA and GGA. For periodic systems, meta-GGA functionals are usable as well. To specify the exchange-correlation potentials of Libxc library, there are two ways. If the exchange and correlation potentials are given separately, you need to specify both ``alibx`` and ``alibc`` separately. If the exchange and correlation potentials are given as a combined set, you need to specify ``alibxc``. We show below an example:
   |    &functional
   |       alibx = 'LDA_X'
   |       alibc = 'LDA_C_PZ'
   | Note that, the hybrid functionals (hybrid gga/mgga) are not supported.
   |
   | To use libxc libraries, ``--enable-libxc`` option must be added in excecuting configure. The available option of the exchange-correlation functionals are listed in the LibXC website. [See http://www.tddft.org/programs/libxc/functionals/]

.. _&rgrid:

&rgrid
------

.. _dl(3):

dl(3)
^^^^^

real(8), default=0d0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Spacing of real-space grids.
   | This cannot be used together with ``&rgrid/num_rgrid``.

.. _num_rgrid(3):

num_rgrid(3)
^^^^^^^^^^^^

integer, default=0

   | Available for the DFT/TDDFT based options of ``theory``.
   | Number of real-space grids for each direction.
   | This cannot be used together with ``&rgrid/dl``.

.. _&kgrid:

&kgrid
------

.. _num_kgrid(3):

num_kgrid(3)
^^^^^^^^^^^^

integer, default=1

   | Available for ``yn_periodic='y'`` in the DFT/TDDFT based options of ``theory``.
   | Number of k-points (grid points of k-vector) for each direction discretizing the Brillouin zone.

.. _file_kw:

file_kw
^^^^^^^

character, default='none'

   | Available for ``yn_periodic='y'`` in the DFT/TDDFT based options of ``theory``.
   | File name for a file that includes user specified k-points. This file will be read if ``num_kgrid`` is equal to 0 or negative values. The file should be described in the following format :
   |
   |   8     #(number of k-points)
   |   1   -0.50  -0.50  -0.50   0.1250   #(id, kx, ky, kz, weight)
   |   2   -0.50  -0.50   0.00   0.1250
   |   3   -0.50   0.00  -0.50   0.1250
   |   4   -0.50   0.00   0.00   0.1250
   |   5    0.00  -0.50  -0.50   0.1250
   |   6    0.00  -0.50   0.00   0.1250
   |   7    0.00   0.00  -0.50   0.1250
   |   8    0.00   0.00   0.00   0.1250

.. /&tgrid:

&tgrid
------

.. _nt:

nt
^^

integer, Default=0

   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   | Number of total time steps for real-time propagation.

.. _dt:

dt
^^

real(8), Default=0d0

   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   | Time step size.

.. _gram_schmidt_interval:

gram_schmidt_interval
^^^^^^^^^^^^^^^^^^^^^

integer, default=-1

   | Available for TDDFT-based options of ``theory``.
   | Interval of a time step for the Gram-Schmidt orthonormalization of the orbitals during time evolution calculations. If this is set to a negative value, no Gram-Schmidt orthogonalization will be achieved. If this is set to zero, the Gram-Schumidt orthogonalization is carried out once at the initial step only. Usually this Gram-Schmidt orthogonalization is not necessary and should not be used.

.. _&propagation:

&propagation
------------

.. _n_hamil:

n_hamil
^^^^^^^

integer, default=4
   | Available for TDDFT-based options of ``theory``.
   | Order of the Taylor expansion adopted for the propagation operator.

.. _propagator:

propagator
^^^^^^^^^^

character, default=middlepoint

   | Available for TDDFT-based options of ``theory``.
   | Choice of the propagator in the time evolution calculation.
   | Options:
   |   ``middlepoint`` / Hamiltoinan at midpoint of two-times is used in the propagation if ``yn_predictor_corrector = 'y'``. Hamiltoian at the time :math:`t` is used if ``yn_predictor_corrector = 'n'``.
   |   ``aetrs`` / time-reversal symmetry propagator. [M.A.L. Marques, A. Castro, G.F. Bertsch, and A. Rubio, Comput. Phys. Commun., 151 60 (2003)].

.. _yn_predictor_corrector:

yn_predictor_corrector
^^^^^^^^^^^^^^^^^^^^^^

character, default='n'
   | Available for TDDFT-based options of ``theory``.
   | Switch of the predictor-corrector method of TDDFT.
   | For meta-GGA functionals (``xc='tbmbj'``), the predictor corrector is automatically used even with ``yn_predictor_corrector='n'``.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_fix_func:

yn_fix_func
^^^^^^^^^^^

character, default='n'
   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   | Switch not to update the Hamiltonian during the time evolution, i.e., ground state Hamiltonian is used during the propagation.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _&scf:

&scf
----

.. _method_init_wf:

method_init_wf
^^^^^^^^^^^^^^

character, default='gauss'

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | The generation method of the initial orbitals at the begening of the SCF iteration in DFT calculations. For a stable calculation of very large systems, multiple gaussian functions are preferred for a stable calculation.
   | Options:
   |   ``gauss`` / single gauss function per orbital centered at a position determined by random numbers
   |   ``gauss2`` / two gauss functions per orbital centered at positions determined by random numbers
   |   ``gauss3`` / three gauss functions per orbital centered at positions determined by random numbers
   |   ``gauss4`` / four gauss functions per orbital centered at positions determined by random numbers
   |   ``gauss5`` / five gauss functions per orbital centered at positions determined by random numbers
   |   ``gauss10`` / ten gauss functions per orbital centered at positions determined by random numbers
   |   ``random`` / a random number is assigned at each real-space grid point of orbitals

.. _method_init_density:

method_init_density
^^^^^^^^^^^^^^^^^^^

[Trial] character, default='wf'

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Specifying how to generate the initial density to start the SCF iteration in the DFT calculation.
   | Supported for limited formats of pseudopotentials ('KY' and 'UPF').
   | Options:
   |  ``wf`` / generate from the initial wavefunctions (cf. ``method_init_wf``).
   |  ``pp`` / generate from a superposition of the pseudo-atom densities.

.. _iseed_number_change:

iseed_number_change
^^^^^^^^^^^^^^^^^^^

integer, default=0

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Change a seed of random numbers that are used to generate initial orbitals. The value specified by this parameter is added to the seed.

.. _nscf:

nscf
^^^^

integer, Default=300

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Number of maximum SCF iterations in the DFT calculation.

.. _method_min:

method_min
^^^^^^^^^^

character, Default='cg'

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Method for updating orbitals in the SCF iteration. At present only confjugate gradient method is implemented.
   | Options:
   |  ``cg`` / Conjugate-Gradient(CG) method

.. _ncg:

ncg
^^^

integer, default=4

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Number of interations of conjugate-gradient method in the SCF iteration.

.. _ncg_init:

ncg_init
^^^^^^^^

integer, default=4

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Number of interations of conjugate-gradient method for the first SCF step.

.. _method_mixing:

method_mixing
^^^^^^^^^^^^^

character, default='broyden'

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Method to update density/potential in the scf iteration.
   | Options:
   |  ``simple`` / Simple mixing method
   |  ``broyden`` / modified Broyden method
   |  ``pulay`` / Pulay method

.. _mixrate:

mixrate
^^^^^^^

real(8), default=0.5d0

   | Available for ``method_mixing='simple'`` in 'dft' and 'dft_md' options of ``theory``.
   | Mixing ratio for simple mixing.

.. _nmemory_mb:

nmemory_mb
^^^^^^^^^^

integer, default=8

   | Available for ``method_mixing='broyden'`` in 'dft' and 'dft_md' options of ``theory``.
   | Number of previous densities to be stored in the SCF iteration using the modified Broyden method. This must be less than 21.

.. _alpha_mb:

alpha_mb
^^^^^^^^

real(8), default=0.75d0

   | Available for ``method_mixing='broyden'`` in 'dft' and 'dft_md' options of ``theory``.
   | A parameter of the modified Broyden method.

.. _nmemory_p:

nmemory_p
^^^^^^^^^

integer, default=4

   | Available for ``method_mixing='pulay'`` in 'dft' and 'dft_md' options of ``theory``.
   | Number of previous densities to be stored in the SCF iteration using the Pulay method.

.. _beta_p:

beta_p
^^^^^^

real(8), default=0.75d0

   | Available for ``method_mixing='pulay'`` in 'dft' and 'dft_md' options of ``theory``.
   | A parameter of the mixing rate of the Pulay method.

.. _yn_auto_mixing:

yn_auto_mixing
^^^^^^^^^^^^^^

character, default='n'

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Switch to change the mixing rate automatically (i.e. automatic adjustments of ``mixrate``/``alpha_mb``/``beta_p``)
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _update_mixing_ratio:

update_mixing_ratio
^^^^^^^^^^^^^^^^^^^

real(8), default=3.0d0

   | Available for ``yn_auto_mixing='y'`` in 'dft' and 'dft_md' options of ``theory``.
   | Threshold for the change of the mixing rate in ``yn_auto_mixing='y'`` option. The mixing-rate is reduced to half when the ratio of the density differences between the current and previous iteration steps is larger than ``update_mixing_ratio``.

.. _yn_subspace_diagonalization:

yn_subspace_diagonalization
^^^^^^^^^^^^^^^^^^^^^^^^^^^

character, default='y'

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Switch for the subspace diagonalization during SCF iterations.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _convergence:

convergence
^^^^^^^^^^^

character, default='rho_dne'

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Specify a quantity that is used for convergence check of the SCF iteration.
   | Options:
   |   ``'rho_dne'``/ Convergence is checked by sum_ix|rho(ix,iter)-rho(ix,iter-1)|dx/N. N is ``&system/nelec``.
   |   ``'norm_rho'``/ Convergence is checked by the square of the norm of the density difference, ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\=sum_ix|rho(ix,iter)-rho(ix,iter-1)|\ :sup:`2`\.
   |   ``'norm_rho_dng'``/ Convergence is checked by ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\/(number of grids). "dng" means "devided by number of grids".
   |   ``'norm_pot'``/ Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\, where Vlocal is Vh + Vxc + Vps_local.
   |   ``'pot_dng'``/ Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\/(number of grids).

.. _threshold:

threshold
^^^^^^^^^

real(8), default=1d-17 [a.u.] (for ``convergence='rho_dne'``) and -1 (for other options of ``convergence``))

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Threshold of convergence that is specified by ``convergence`` keyword.

.. _nscf_init_redistribution:

nscf_init_redistribution
^^^^^^^^^^^^^^^^^^^^^^^^

integer, default=10

   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Number of initial iterations during which a redistribution of the occupation number is suppressed in the finite temperature calculation.

.. _nscf_init_no_diagonal:

nscf_init_no_diagonal
^^^^^^^^^^^^^^^^^^^^^

integer, default=10

   | Available for ``&scf/yn_subspace_diagonalization='y'`` in 'dft' option of ``theory``.
   | Number of initial iterations during which the subspace diagonalization will not be carried out.

.. _nscf_init_mix_zero:

nscf_init_mix_zero
^^^^^^^^^^^^^^^^^^

integer, default=-1

   | Available for 'dft' option of ``theory``.
   | The density will not be mixed (i.e. fixed) during the given number of the SCF iteration, that is, orbitals are optimized without updating the density.

.. _conv_gap_mix_zero:

conv_gap_mix_zero
^^^^^^^^^^^^^^^^^

real(8), default=99999d0

   | Available if ``nscf_init_mix_zero`` is positive value in the 'dft' option of ``theory``.
   | Specify a condition to quit the fixed density iteration forced by ``step_initial_mix_zero`` option. Mixing of the density will start after the band-gap energy exceeds this parameter for consecutive five SCF iteration steps.

.. _&emfield:

&emfield
--------

.. _trans_longi:

trans_longi
^^^^^^^^^^^

character, default='tr'

   | Available for ``yn_periodic='y'`` in 'maxwell' and TDDFT based options of ``theory``.
   | Specify the treatment of the polarization in the time evolution calculation.
   | Options:
   |   ``'tr'`` / Transverse
   |   ``'lo'`` / longitudinal
   |   ``'2d'`` / 2D maxwell-TDDFT (2D approximation) method (for more details, see ``film_thickness`` of &maxwell)

.. _ae_shape1:

ae_shape1
^^^^^^^^^

.. _ae_shape2:

ae_shape2
^^^^^^^^^

character, Default='none'

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Envelope shape of the first/second pulse. 'Acos2' indicates a cosine square envelope for vector potential, and 'Ecos2' a cosine square envelope for electric field.
   | Options:
   |   ``'impulse'`` / A weak impulsive field is applied at :math:`t=0`. This will be used to explore linear response properties. The magnitude of the impulse can be specified by ``e_impulse``.
   |   ``'Acos2'`` / Envelope of cos\ :sup:`2`\ for a vector potential.
   |   ``'Acos3'`` / Envelope of cos\ :sup:`3`\ for a vector potential.
   |   ``'Acos4'`` / Envelope of cos\ :sup:`4`\ for a vector potential.
   |   ``'Acos6'`` / Envelope of cos\ :sup:`6`\ for a vector potential.
   |   ``'Acos8'`` / Envelope of cos\ :sup:`8`\ for a vector potential.
   |   ``'Ecos2'`` / Envelope of cos\ :sup:`2`\ for an electric field.
   |   ``'Asin2cos'`` [Trial] / Envelope of sin\ :sup:`2`\ with cosine type oscillation for a vector potential.
   |   ``'Asin2_cw'`` [Trial] / Envelope of sin\ :sup:`2`\ at the beginning and continuous wave after that for a vector potential (for 'ae_shape1' only).
   |   ``'input'`` [Trial] / read the vector potential as a numerical table with ``file_input1`` option (for 'ae_shape1' only).
   |   ``'none'`` / no incident field is applied.
   |
   | If 'Ecos2' is adopted, 'phi_cep1' must be chosen either 0.75 or 0.25, since otherwise the time integral of the electric field (vector potential at the end of the pulse) does not vanishi. There is no such restriction for 'Acos2' pulses.
   |
   | For ``yn_periodic='n'``, available choices are limited to ``'impulse'``, ``'Acos2'``, and ``'Ecos2'``.

.. _file_input1:

file_input1
^^^^^^^^^^^

character, default=''

   | Available if ``ae_shape1='input'`` is specified and ``theory='tddft_pulse'``.
   | Name of an input file that contains user-defined vector potential. The file must be a numerical table separated by blank, having four columns; the first column is time and second to fourth columns are Ax/c, Ay/c, Az/c, repsectively. All the quantities are written using the units specified by ``unit_system``. '#' and '!' may be used for a comment line.
   | Note that a linear interpolation will be applied when the time step differs from that used in the calculation.

.. _e_impulse:

e_impulse
^^^^^^^^^

real(8), Default=1d-2 a.u.

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Magnitude of the impulse in the impulsive perturbation. This valiable has the dimention of momentum, energy*time/length.

..
  #(commented out: not implemented yet)
  #- **t_impulse**
  #   | Available for ``theory='XXX'``.
  #   not yet implemented XXX
..

.. _E_amplitude1:

E_amplitude1
^^^^^^^^^^^^

.. _E_amplitude2:

E_amplitude2
^^^^^^^^^^^^

real(8), default=0d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Maximum amplitude of electric field for the first/second pulse. This valiable has the dimension of electric field, energy/(length*charge). This cannot be set with ``&emfield/I_wcm2_1`` (``I_wcm2_2``) simultaneously.

.. _I_wcm2_1:

I_wcm2_1
^^^^^^^^

.. _I_wcm2_2:

I_wcm2_2
^^^^^^^^

real(8), default=-1d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Maximum intensity (W/cm\ :sup:`2`\) of the first/second pulse. This valiable cannot be set with ``&emfield/E_amplitude1`` (``E_amplitude2``) simultaneously. For this quantity, a unit of W/cm\ :sup:`2`\ is adopted irrespective of ``&units\unit_system``.

.. _tw1:

tw1
^^^

.. _tw2:

tw2
^^^

real(8), default=0d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Duration of the first/second pulse (edge-to-edge time length).
   | Note that this is not the FWHM duration.

.. _omega1:

omega1
^^^^^^

.. _omega2:

omega2
^^^^^^

real(8), default=0d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Mean photon energy (average frequency multiplied by the Planck constant) of the first/second pulse.

.. _epdir_re1(3):

epdir_re1(3)
^^^^^^^^^^^^

.. _epdir_re2(3):

epdir_re2(3)
^^^^^^^^^^^^

real(8), default=1d0, 0d0, 0d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Real part of the polarization unit vector for the first/second pulse.

.. _epdir_im1(3):

epdir_im1(3)
^^^^^^^^^^^^

.. _epdir_im2(3):

epdir_im2(3)
^^^^^^^^^^^^

real(8), default=0d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Imaginary part of the polarization unit vector for the first/second pulse. Using both real 'epdir_re1' and imaginary 'epdir_im1' parts of the polarization vector, circularly and general ellipsoidary polarized pulses may be described.

.. _phi_cep1:

phi_cep1
^^^^^^^^

.. _phi_cep2:

phi_cep2
^^^^^^^^

real(8), default=0d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Carrier envelope phase of the first/second pulse. It specifies the CEP in unit of :math:`2\pi`.

.. _t1_t2:

t1_t2
^^^^^

real(8), default=0d0
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Time-delay between the first and the second pulses.

.. _t1_start:

t1_start
^^^^^^^^

real(8), default=0d0

   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Shift the starting time of the first pulse. (this is not available for multiscale option).

.. _num_dipole_source:

num_dipole_source
^^^^^^^^^^^^^^^^^

integer, default=0

   | Available for TDDFT based options of ``theory``.
   | Number of radiation sources to mimic optical near fields. Maximum number is ``2``.

.. _vec_dipole_source(3,num_dipole_source):

vec_dipole_source(3,num_dipole_source)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for TDDFT based options of ``theory``.
   | Dipole vectors of the radiation sources mimicing optical near fields.

.. _cood_dipole_source(3,num_dipole_source):

cood_dipole_source(3,num_dipole_source)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for TDDFT based options of ``theory``.
   | Coordinates of the radiation sources mimicing optical near fields.

.. _rad_dipole_diele:

rad_dipole_diele
^^^^^^^^^^^^^^^^

real(8), default=2d0 [a.u.]

   | Available for TDDFT based options of ``theory``.
   | Radii of dielectric spheres of the radiation sources mimicing optical near fields.

.. _&singlescale[Trial]:

&singlescale[Trial]
-------------------

.. _method_singlescale:

method_singlescale
^^^^^^^^^^^^^^^^^^

character, default='3d'

   | Available for ``theory='single_scale_maxwell_tddft'``.
   | Type of single-scale Maxwell-TDDFT method.
   | Options:
   | ``'3d'`` / 3-dimensional FDTD + TDDFT
   | ``'1d'`` / 1-dimensional FDTD (along the z axis) + TDDFT
   | ``'1d_fourier'`` / ``'1d'`` with 3D Fourier component of the vector potential

.. _cutoff_G2_emfield:

cutoff_G2_emfield
^^^^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``theory='single_scale_maxwell_tddft'``.
   | Cutoff energy of Fourier component of the vector potential when method_singlescale='1d_fourier'.

.. _yn_symmetrized_stencil:

yn_symmetrized_stencil
^^^^^^^^^^^^^^^^^^^^^^

[Trial] character, default='n'

   | Available for ``theory='single_scale_maxwell_tddft'``.
   | Switch to symmetrize the finite differences for the product of vector potential and orbitals, :math:`(\nabla A(r) \cdot \psi(r))`. This option improves hermiticity of the Hamiltonian although computational cost increases.

.. _yn_put_wall_z_boundary:

yn_put_wall_z_boundary
^^^^^^^^^^^^^^^^^^^^^^

[Trial] character, default='n'

   | Available for DFT/TDDFT based options of ``theory``.
   | Option to put potential wall near the boundary planes at *z*\ =0 and *z*\ =``&system/al(3)``. This potential prevents electrons from crossing the *z*\ -boundary plane. In the single-scale Maxwell-TDDFT method, the electron density on the *z*\ -boundary plane harms the norm conservation of electrons due to the discontinuity of the vectorpotential. The wall is described using the square of cosine function.
   | Options:
   |   ``'y'`` / put the potential wall
   |   ``'n'`` / no potential wall

.. _wall_height:

wall_height
^^^^^^^^^^^

real(8), default=100.0 [eV]

   | Available for ``yn_put_wall_z_boundary='y'``.
   | The height of the potential wall.

.. _wall_width:

wall_width
^^^^^^^^^^

real(8), default=5.0 [Angstrom]

   | Available for ``yn_put_wall_z_boundary='y'``.
   | The width of the potential wall defined by the length from the potential peak (\ *z*\ =0 and *z*\ =``&system/al(3)``) to the edge.

.. _&multiscale:

&multiscale
-----------

.. _fdtddim:

fdtddim
^^^^^^^

[Trial] character, default='1d'

   | Available for ``theory='multi_scale_maxwell_tddft'`` with ``yn_periodic='y'`` and ``theory='maxwell_sbe'``
   | Dimension of macroscopic scale system (Maxwell(FDTD) calculation) in multiscale Maxwell-TDDFT method.
   | Options:
   | ``'3d'`` / 3-dimensional FDTD for macroscopic electromagnetism [currently not available]
   | ``'1d'`` / 1-dimensional FDTD (along the *x*\ -axis) for macroscopic electromagnetism

.. _nx_m:

nx_m
^^^^

integer, default=1

   | Available for ``theory='multi_scale_maxwell_tddft'`` with ``yn_periodic='y'`` and ``theory='maxwell_sbe'``
   | Number of macroscopic grid points inside materials for *x*\ -direction.

.. _ny_m:

ny_m
^^^^

integer, default=1

   | Available for ``theory='multi_scale_maxwell_tddft'`` with ``yn_periodic='y'`` and ``theory='maxwell_sbe'``
   | Number of macroscopic grid points inside materials for *y*\ -direction.

.. _nz_m:

nz_m
^^^^

[Trial] integer, default=1)

   | Available for ``theory='multi_scale_maxwell_tddft'`` with ``yn_periodic='y'`` or ``theory='maxwell_sbe'``
   | Number of macroscopic grid points inside materials for (\ *y*\ /\ *z*\ )-direction.

.. _hx_m:

hx_m
^^^^

real(8), default=0d0
   | Available for ``theory='multi_scale_maxwell_tddft'`` with ``yn_periodic='y'`` or ``theory='maxwell_sbe'``
   | Grid spacing of macroscopic coordinate for *x*\ -direction.
   | Variable ``hx_m`` is deprecated, and will be moved to ``&units/dl_em(1)``

.. _hy_m:

hy_m

.. _hz_m:

hz_m
^^^^

[Trial] real(8), default=0d0

   | Available for ``theory='multi_scale_maxwell_tddft'`` with ``yn_periodic='y'`` or ``theory='maxwell_sbe'``
   | Grid spacing of macroscopic coordinate for (\ *y*\ /\ *z*\ )-direction.
   | Variable ``hy_m`` and ``hz_m`` are deprecated, and will be moved to ``&units/dl_em(2:3)``

.. _nxvacl_m:

nxvacl_m
^^^^^^^^
integer, default=1/0

   | Available for ``theory='multi_scale_maxwell_tddft'`` or ``'maxwell_sbe'``
   | The parameter ``nxvacl_m`` will be replaced by ``nxvac_m`` and eventually removed in the future.

.. _nxvacr_m:

nxvacr_m
^^^^^^^^

integer, default=1/0

   | Available for ``theory='multi_scale_maxwell_tddft'`` or ``'maxwell_sbe'``
   | The parameter ``nxvacr_m`` will be replaced by ``nxvac_m`` and eventually removed in the future.

.. _nxvac_m(2):

nxvac_m(2)
^^^^^^^^
integer, default=0

   | Available for ``theory='multi_scale_maxwell_tddft'`` or ``'maxwell_sbe'``
   | Represents the number of vacuum cells between the edge of the material region and the computational boundary. The first element of the array represents the number of cells from the leftmost cell (ix=1) on the x-axis to the left boundary. The second element represents the number of cells from the rightmost cell on the x-axis to the right boundary.

nyvac_m(2)
^^^^^^^^
integer, default=0

   | Available for ``theory='multi_scale_maxwell_tddft'`` or ``'maxwell_sbe'``
   | Provides same functionality of ``nxvac_m(2)`` for y-direction.

nzvac_m(2)
^^^^^^^^
integer, default=0

   | Available for ``theory='multi_scale_maxwell_tddft'`` or ``'maxwell_sbe'``
   | Provides same functionality of ``nxvac_m(2)`` for z-direction.


.. _&maxwell:


&maxwell
--------

.. _al_em(3):

al_em(3)
^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | Size of simulation box in electromagnetic analysis.
   | Only two of ``al_em``, ``dl_em``, and ``num_rgrid_em`` must be set.

.. _dl_em(3):

dl_em(3)
^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'`` and  ``theory='multi_scale_maxwell_tddft'``.
   | Spacing of real-space grids in electromagnetic analysis.
   | Only two of ``al_em``, ``dl_em``, and ``num_rgrid_em`` must be set.

.. _num_rgrid_em(3):

num_rgrid_em(3)
^^^^^^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'``.
   | Number of real-space grids in electromagnetic analysis.
   | Only two of ``al_em``, ``dl_em``, and ``num_rgrid_em`` must be set.

.. _at_em:

at_em
^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | Total time for electromagnetic analysis.
   | Two of ``at_em``, ``dt_em``, and ``nt_em`` must be set.
   | Otherwise, both ``at_em`` and ``nt_em`` or either of those must be set.
   | (For the latter, ``dt_em`` is automatically determined from CFL condition)

.. _dt_em:

dt_em
^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | Time step size for electromagnetic analysis.
   | If default is selected, this is automatically determined from CFL condition.

.. _nt_em:

nt_em
^^^^^

integer, default=0

   | Available for ``theory='maxwell'``.
   | Number of total time steps of time propagation in electromagnetic analysis.

.. _boundary_em(3,2):

boundary_em(3,2)
^^^^^^^^^^^^^^^^

character, default='default'

   | Available for ``theory='maxwell'`` and ``theory='multi_scale_maxwell_tddft'`` and ``theory='maxwell_sbe'``
   | Boundary condition in electromagnetic analysis. The first index(1-3 rows) corresponds to *x*\ , *y*\ , and *z* axes. The second index(1-2 columns) corresponds to bottom and top of the axes.
   | Options:
   | ``'abc'`` / absorbing boundary
   | ``'pec'`` / perfect electric conductor
   | ``'periodic'`` / periodic boundary
   |
   | If ``&system/yn_periodic='n'``, ``'default'``, ``'abc'``, and ``'pec'`` can be chosen, where ``'default'`` automatically chooses ``'abc'``. If ``&system/yn_periodic='y'``, ``'default'``, ``'abc'``, and ``'periodic'`` can be chosen, where ``'default'`` automatically chooses ``'periodic'``. | When ``theory='maxwell'``, perfectly matched layer(PML) is used for ``'abc'``.

.. _shape_file:

shape_file
^^^^^^^^^^

character, default='none'

   | Available for ``theory='maxwell'``.
   | Name of input shape file in electromagnetic analysis. The shape file can be generated by using ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).

.. _media_num:

media_num
^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'``  and ``theory='maxwell_sbe'``.
   | Number of media in electromagnetic analysis.

.. _media_type(:):

media_type(:)
^^^^^^^^^^^^^

character, default='vacuum'

   | Available for ``theory='maxwell'`` and ``theory='maxwell_sbe'``
   | ``media_type(n)`` spesifies type of n-th media in electromagnetic analysis.
   | Options:
   |   ``'vacuum'``
   |   ``'constant media'``
   |   ``'pec'``
   |   ``'lorentz-drude'``
   | If ``'lorentz-drude'`` is chosen, linear response calculation is feasible by setting ``&emfield/ae_shape1 or ae_shape2='impulse'``.
   |  Besides, in the case of  ``theory='maxwell_sbe'``, ``'multiscale'`` is also available.

.. _epsilon_em(:):

epsilon_em(:)
^^^^^^^^^^^^^

real(8), Default=1d0

   | Available for ``theory='maxwell'``, ``theory='maxwell_sbe'`` and for TDDFT based options of ``theory`` with ``trans_longi='2d'``.
   | For ``theory='maxwell'``, ``epsilon_em(n)`` spesifies relative permittivity of n-th media in electromagnetic analysis.
   | For TDDFT based options of ``theory`` with ``trans_longi='2d'``, the relative permittivity of the transparent media on both sides of the film is specified by ``epsilon_em(1)`` and ``epsilon_em(2)``, respectively.

.. _mu_em(:):

mu_em(:)
^^^^^^^^

real(8), default=1d0

   | Available for ``theory='maxwell'``.
   | ``mu_em(n)`` spesifies relative permeability of n-th media in electromagnetic analysis.

.. _sigma_em(:):

sigma_em(:)
^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | ``sigma_em(n)`` spesifies conductivity of n-th media in electromagnetic analysis.

.. _pole_num_ld(:):

pole_num_ld(:)
^^^^^^^^^^^^^^

integer, default=1

   | Available for ``theory='maxwell'``.
   | ``pole_num_ld(n)`` spesifies number of poles of n-th media, available for ``type_media(n)='lorentz-drude'`` in electromagnetic analysis.

.. _omega_p_ld(:):

omega_p_ld(:)
^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | ``omega_p_ld(n)`` spesifies plasma frequency of n-th media, available for ``type_media(n)='lorentz-drude'`` in electromagnetic analysis.

.. _f_ld(:,:):

f_ld(:,:)
^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | ``f_ld(n,m)`` spesifies m-th oscillator strength of n-th media, available for ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is the media ID whose maximum value is given by ``media_num``. The second index is the pole ID whose maximum value is given by ``pole_num_ld(n)``.

.. _gamma_ld(:,:):

gamma_ld(:,:)
^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | ``gamma_ld(n,m)`` spesifies m-th collision frequency of n-th media, available for ``type_media(n)='lorentz-drude'`` in electromagnetic analysis. The first index is the media ID whose maximum value is given by ``media_num``. The second index is the pole ID whose maximum value is given by ``pole_num_ld(n)``.

.. _omega_ld(:,:):

omega_ld(:,:)
^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | ``omega_ld(n,m)`` spesifies m-th oscillator frequency of n-th media, available for ``type_media(n)='lorentz-drude'`` in electromagnetic analysis. The first index is the media ID whose maximum value is given by ``media_num``. The second index is the pole ID whose maximum value is given by ``pole_num_ld(n)``.

.. _wave_input:

wave_input
^^^^^^^^^^

character, default='none'

   | Available for ``theory='maxwell'``.
   | If ``'source'``, the incident pulse in electromagnetic analysis is generated by the incident current source.

.. _ek_dir1(3):

ek_dir1(3)
^^^^^^^^^^

.. _ek_dir2(3):

ek_dir2(3)
^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | Propagation direction of the first/second pulse (\ *x*\ , *y*\ , and *z* directions). Each component must be 0d0 or 1d0.

.. _source_loc1(3):

source_loc1(3)
^^^^^^^^^^^^^^

.. _source_loc2(3):

source_loc2(3)
^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | Location of the incident current source of the first/second pulse. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/yn_periodic='n'`` while ranges from ``0`` to ``al_em`` for ``&system/yn_periodic='y'``.

.. _gbeam_sigma_plane1(3):

gbeam_sigma_plane1(3)
^^^^^^^^^^^^^^^^^^^^^

.. _gbeam_sigma_plane2(3):

gbeam_sigma_plane2(3)
^^^^^^^^^^^^^^^^^^^^^

.. _gbeam_sigma_line1(3):

gbeam_sigma_line1(3)
^^^^^^^^^^^^^^^^^^^^
.. _gbeam_sigma_line2(3):

gbeam_sigma_line2(3)
^^^^^^^^^^^^^^^^^^^^

[Trial] real(8), default=-1d0

   | Available for ``theory='maxwell'`` with ``wave_input='source'``.
   | These input keywords specify the width of Gauss function, exp(-0.5(abs(r-r_0)/sigma)^2), applied for the incident current source to generate the first/second pulse. These input keywords work only when their values > 0.0d0. The center of the Gauss function, r_0, is specified by ``source_loc1/2``. ``gbeam_sigma_plane1/2`` specifies the width of 2D Gauss function (\ *xy*\ , *yz*\ , and *xz* planes). ``gbeam_sigma_line1/2`` specifies the width of 1D Gauss function (\ *x*\ , *y*\ , and *z* axes).

.. _obs_num_em:

obs_num_em
^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'``.
   | Number of observation points in electromagnetic analysis. From the obtained results, figure and animation files can be generated by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

.. _obs_samp_em:

obs_samp_em
^^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'``.
   | Sampling time-step of the observation in electromagnetic analysis.
   | If default is selected, this is automatically determined.

.. _obs_loc_em(:,3):

obs_loc_em(:,3)
^^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'``.
   | ``obs_loc_em(n,1:3)=x,y,z`` spesifies location of the n-th observation point in electromagnetic analysis. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/yn_periodic='n'`` while ranges from ``0`` to ``al_em`` for ``&system/yn_periodic='y'``.

.. _obs_plane_ene_em(:,:):

obs_plane_ene_em(:,:)
^^^^^^^^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``theory='maxwell'``.
   | ``obs_loc_em(n,:)=energy1,energy2,energy3,...`` spesifies energy value of the n-th observation point in electromagnetic analysis. At the spesified energies, Fourier-transformed spatial distributions on the xy, yz, and xz plans are outputed. This input keyword must be larger than 0.

.. _yn_obs_plane_em(:):

yn_obs_plane_em(:)
^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='maxwell'``.
   | Spesify whether or not to generate output of the electrmagnetic fields on the planes (\ *xy*\ , *yz*\ , and *xz* planes) for n-th observation point. This option must be ``'y'`` for generating animation files by using ``FDTD_make_figani`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).
   | Options:
   |   ``'y'``
   |   ``'n'``

.. _yn_obs_plane_integral_em(:):

yn_obs_plane_integral_em(:)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='maxwell'``.
   | Specify whether or not to generate output of the spatial integration of electromagnetic fields on the planes (\ *xy*\ , *yz*\ , and *xz* planes) for n-th observation point.
   | Options:
   |   ``'y'``
   |   ``'n'``

.. _yn_wf_em:

yn_wf_em
^^^^^^^^

character, default='y'

   | Available for ``theory='maxwell'``.
   | Switch of a window function for linear response calculation.
   | Options:
   |   ``'y'``
   |   ``'n'``

.. _film_thickness:

film_thickness
^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for TDDFT based options of ``theory`` with ``trans_longi='2d'``.
   | Thickness of the film for the 2D maxwell-TDDFT (2D approximation) method [S. Yamada and K. Yabana, PRB 103, 155426 (2021)].
   | For a slab system, ``film_thickness`` should be set to the side length of the calculation cell, i.e., the thickness of the slab plus the length of the vacuum region [S. Yamada et al., PRB 98, 245147 (2018)].
   | The relative permittivity of the transparent media on both sides of the film can be specified by ``epsilon_em(1)`` and ``epsilon_em(2)``, respectively (default=vacuum).

.. _media_id_pml(3:2):

media_id_pml(3:2)
^^^^^^^^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'``.
   | Media ID used in PML. The first index(1-3 rows) corresponds to *x*\ , *y*\ , and *z* axes. The second index(1-2 columns) corresponds to bottom and top of the axes.

.. _media_id_source1:

media_id_source1
^^^^^^^^^^^^^^^^

.. _media_id_source2:

media_id_source2
^^^^^^^^^^^^^^^^

integer, default=0
   | Available for ``theory='maxwell'``.
   | Media ID used in incident current source1/source2 to generate the first/second pulse.

.. _bloch_k_em:

bloch_k_em(3)
^^^^^^^^^^^^^

[Trial] real(8), default=0d0

   | Available for ``theory='maxwell'`` with ``yn_periodic='y'``.
   | Wavenumber used in Bloch boundary conditions. When sum(|bloch_k_em(:)|)>0, Bloch boundary conditions are automatically applied.

.. _bloch_k_em:

bloch_real_imag_em(3)
^^^^^^^^^^^^^^^^^^^^^

[Trial] character, default='real'

   | Available for ``theory='maxwell'`` with ``yn_periodic='y'`` and sum(|bloch_k_em(:)|)>0.
   | Specify real or imaginary parts for exp(ikr) used in Bloch boundary conditions.
   | Options:
   |   ``'real'``
   |   ``'imag'``

.. _ase_num_em:

ase_num_em
^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'`` with ``yn_periodic='n'``.
   | Number of energy or wavelength grid points specified by ``ase_ene_min_em/ase_ene_max_em`` or ``ase_wav_min_em/ase_wav_max_em``.
   | If this is specified as larger than 0, Absorption-, Scattering-, and Extinction-cross-sections will be outputed at the end of calculation.
   | Those are normalized by the spectral distribution of the incident pulse.

.. _ase_ene_min_em:

ase_ene_min_em
^^^^^^^^^^^^^^
.. _ase_ene_max_em:

ase_ene_max_em
^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``theory='maxwell'`` with ``ase_num_em>0`` and ``yn_periodic='n'``.
   | Energy range for Absorption-, Scattering-, and Extinction-cross-sections.

.. _ase_wav_min_em:

ase_wav_min_em
^^^^^^^^^^^^^^
.. _ase_wav_max_em:

ase_wav_max_em
^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``theory='maxwell'`` with ``ase_num_em>0`` and ``yn_periodic='n'``.
   | Wavelength range for Absorption-, Scattering-, and Extinction-cross-sections.

.. _ase_smedia_id_em:

ase_smedia_id_em
^^^^^^^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'`` with ``ase_num_em>0`` and ``yn_periodic='n'``.
   | Media ID used as surrounding media.

.. _ase_box_cent_em(3):

ase_box_cent_em(3)
^^^^^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'`` with ``ase_num_em>0`` and ``yn_periodic='n'``.
   | ``ase_box_cent_em(1:3)=x,y,z`` spesifies location of the center of a closed surface (box shape) to calculate Absorption-, Scattering-, and Extinction-cross-sections.

.. _ase_box_size_em(3):

ase_box_size_em(3)
^^^^^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``theory='maxwell'`` with ``ase_num_em>0`` and ``yn_periodic='n'``.
   | ``ase_box_size_em(1:3)=X,Y,Z`` spesifies size of a closed surface (box shape) to calculate Absorption-, Scattering-, and Extinction-cross-sections.

.. _art_num_em:

art_num_em
^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'`` with ``yn_periodic='y'``.
   | Number of energy or wavelength grid points specified by ``art_ene_min_em/art_ene_max_em`` or ``art_wav_min_em/art_wav_max_em``.
   | If this is specified as larger than 0, Absorption-, Reflection-, and Transmission-ratas will be outputed at the end of calculation.
   | Those are normalized by the spectral distribution of the incident pulse.

.. _art_ene_min_em:

art_ene_min_em
^^^^^^^^^^^^^^
.. _art_ene_max_em:

art_ene_max_em
^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``theory='maxwell'`` with ``art_num_em>0`` and ``yn_periodic='y'``.
   | Energy range for Absorption-, Reflection-, and Transmission-ratas.

.. _art_wav_min_em:

art_wav_min_em
^^^^^^^^^^^^^^
.. _art_wav_max_em:

art_wav_max_em
^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``theory='maxwell'`` with ``art_num_em>0`` and ``yn_periodic='y'``.
   | Wavelength range for Absorption-, Reflection-, and Transmission-ratas.

.. _art_smedia_id_em:

art_smedia_id_em
^^^^^^^^^^^^^^^^

integer, default=0

   | Available for ``theory='maxwell'`` with ``art_num_em>0`` and ``yn_periodic='y'``.
   | Media ID used as surrounding media.
   
.. _art_plane_bot_em(3):

art_plane_bot_em(3)
^^^^^^^^^^^^^^^^^^
.. _art_plane_top_em(3):

art_plane_top_em(3)
^^^^^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``theory='maxwell'`` with ``art_num_em>0`` and ``yn_periodic='y'``.
   | ``art_plane_bot_em(1:3)=x1,y1,z1`` and ``art_plane_top_em(1:3)=x2,y2,z2`` spesify location of bottom and top planes on the propagation axis to calculate Absorption-, Reflection-, and Transmission-ratas.

.. _yn_make_shape:

yn_make_shape
^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='maxwell'``.
   | Switch for making shape. This is same functionality for ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).
   | Options:
   |   ``'y'``
   |   ``'n'``

.. _yn_output_shape:

yn_output_shape
^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='maxwell'``.
   | Switch for outputing shape file in cube format when ``yn_make_shape='y'``.
   | Options:
   |   ``'y'``
   |   ``'n'``

.. _yn_copy_x:

yn_copy_x
^^^^^^^^^
.. _yn_copy_y:

yn_copy_y
^^^^^^^^^
.. _yn_copy_z:

yn_copy_z
^^^^^^^^^
character, default='n'

   | Available for ``theory='maxwell'``.
   | See ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).
   | Options:
   |   ``'y'``
   |   ``'n'``

.. _rot_type:

rot_type
^^^^^^^^
character, default='radian'

   | Available for ``theory='maxwell'``.
   | See ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).
   | Options:
   |   ``'radian'``
   |   ``'degree'``

.. _n_s:

n_s
^^^
integer, default=0

   | Available for ``theory='maxwell'``  and ``theory='maxwell_sbe'``.
   | See ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).

.. _typ_s(:):

typ_s(:)
^^^^^^^^
character, default='none'

   | Available for ``theory='maxwell'``  and ``theory='maxwell_sbe'``.
   | See ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).

.. _id_s(:):

id_s(:)
^^^^^^^
integer, default=0

   | Available for ``theory='maxwell'`` and ``theory='maxwell_sbe'``.
   | See ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).

.. _inf_s(:,10):

inf_s(:,10)
^^^^^^^^^^^
real(8), default=0

   | Available for ``theory='maxwell'`` and ``theory='maxwell_sbe'``.
   | See ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).

.. _ori_s(:,3):

ori_s(:,3)
^^^^^^^^^^
.. _rot_s(:,3):

rot_s(:,3)
^^^^^^^^^^
real(8), default=0d0

   | Available for ``theory='maxwell'``
   | See ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html).

.. _&analysis:

&analysis
---------

.. _projection_option:

projection_option
^^^^^^^^^^^^^^^^^

character, default='no'

   | Available for TDDFT based options of ``theory``.
   | Methods of projection to analyze the excited states (e.g. the number of excited electrons).
   | Output files: SYSname_nex.data, SYSname_ovlp.data
   | Options:
   |   ``'no'`` / no projection.
   |   ``'gs'`` / projection to eigenstates of the ground-state Hamiltonian whose k-point is shifted as k+A(t)/c (i.e. Houston functions).
   |   ``'td'`` / projection to instantaneous eigenstates of the time-dependent Hamiltonian.

.. _out_projection_step:

out_projection_step
^^^^^^^^^^^^^^^^^^^

integer, default=100

   | Available when ``projection_option`` is specified.
   | Resuts of the projection analysis will be outputted every ``out_projection_step`` steps during the time-propagation.

.. _threshold_projection:

threshold_projection
^^^^^^^^^^^^^^^^^^^^

real(8), default=1e-6

   | Available when ``projection_option`` is specified.
   | Convergence threshold for the iteration of the eigenstates calculation.z

.. _yn_out_intraband_current:

yn_out_intraband_current
^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available when ``projection_option`` is specified.
   | Switch for output of the intra-band current density [T. Otobe, Phys. Rev. B 94, 235152 (2016).].
   | Output file: SYSname_intra_current.data
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _nenergy:

nenergy
^^^^^^^

integer, default=1000

   | Number of energy grid points for frequency-domain analysis. This parameter is used, for examples, in ``theory='tddft_response'`` and ``theory='maxwell'``.

.. _de:

de
^^

real(8), Default=0.01d0 (eV)

   | Energy grid size for frequency-domain analysis.
   | This parameter is used, for examples, in ``theory='tddft_response'`` and ``theory='maxwell'``.

.. _out_rt_energy_step:

out_rt_energy_step
^^^^^^^^^^^^^^^^^^

integer, default=10

   | Available for the TDDFT based option of ``theory``.
   | Total energy is calculated and printed every ``out_rt_energy_step`` time steps.

.. _yn_out_psi:

yn_out_psi
^^^^^^^^^^

character, default='n'

   | Available for ``theory='dft'``.
   | Switch for output of orbitals.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable
   | The format of the output is specified by &analysis/``format_voxel_data``.

.. _yn_out_dos:

yn_out_dos
^^^^^^^^^^

character, default='n'

   | Available for ``theory='dft'``.
   | Switch for output of density of states.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_out_pdos:

yn_out_pdos
^^^^^^^^^^^

character, default='n'

   | Available for ``theory='dft'``.
   | Switch for output of projected density of states.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_out_dos_set_fe_origin:

yn_out_dos_set_fe_origin
^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available when ``yn_out_dos='y'`` or ``yn_out_pdos='y'``.
   | Switch to set the Fermi energy to zero in the energy axis of DoS.
   | If the temperature is not specified, this option sets the valence band maximum to zero.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_dos_start:

out_dos_start
^^^^^^^^^^^^^

real(8), default=-1d10 (eV)

.. _out_dos_end:

out_dos_end
^^^^^^^^^^^

real(8), default=1d10 (eV)

   | Available when ``yn_out_dos='y'`` or ``yn_out_pdos='y'``.
   | Lower/Upper bound of the energy range for the density of states spectra.
   | If this value is lower/higher than a specific value near the lowest/highest energy level, this parameter is re-set to the value.

.. _out_dos_nenergy:

out_dos_nenergy
^^^^^^^^^^^^^^^

integer, default=601

   | Available when ``yn_out_dos='y'`` or ``yn_out_pdos='y'``.
   | Number of energy points sampled in the density of states spectra.

.. _out_dos_function:

out_dos_function
^^^^^^^^^^^^^^^^

character, default='gaussian'

   | Available when ``yn_out_dos='y'`` or ``yn_out_pdos='y'``.
   | Choice of the smearing function for the density of states spectra.
   | Options:
   |   ``gaussian``  / Gaussian function
   |   ``lorentzian`` / Lorentzian function

.. _out_dos_width:

out_dos_width
^^^^^^^^^^^^^

real(8), default=0.1d0 [eV]

   | Available when ``yn_out_dos='y'`` or ``yn_out_pdos='y'``.
   | Smearing width used in the density of states spectra.

.. _yn_out_dns:

yn_out_dns
^^^^^^^^^^

character, default='n'

   | Available for ``theory='dft'``.
   | Switch to output electron density distribution of the ground state.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_out_dns_rt:

yn_out_dns_rt
^^^^^^^^^^^^^

character, default='n'

   | Available when ``theory='dft_md' or 'theory=tddft_pulse'``.
   | Switch to output electron density distribution during the time-propagation.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_dns_rt_step:

out_dns_rt_step
^^^^^^^^^^^^^^^^^^

integer, default=50

   | Available when ``theory='dft_md' or 'theory=tddft_pulse'``.
   | Density is outputted every ``out_dns_rt_step`` steps.

.. _yn_out_dns_ac_je:

yn_out_dns_ac_je
^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='single_scale_maxwell_tddft'``.
   | Switch to print the electron density, vector potential, electronic current, and ionic coordinates every ``out_dns_ac_je_step`` time steps.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable
   | The data written in binary format are divided into files corresponding to the space-grid parallelization number.

.. _out_dns_ac_je_step:

out_dns_ac_je_step
^^^^^^^^^^^^^^^^^^

integer, default=50

   | Available for ``theory='single_scale_maxwell_tddft'``.
   | Electron density, vector potential, electronic current, and ionic coordinates are outputted every ``outdns_dns_ac_je_step`` time steps.

.. _yn_out_micro_je:

yn_out_micro_je
^^^^^^^^^^^^^^^

character, default='n'

   | Available for TDDFT based methods.
   | Switch to print the microscopic electron current density (``je_micro_***`` files) at every ``out_micro_je_step`` time steps.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_micro_je_step:

out_micro_je_step
^^^^^^^^^^^^^^^^^

integer, default=50

   | Available for TDDFT based methods with ``yn_out_micro_je='y'``.
   | See ``yn_out_micro_je``.


.. _yn_out_dns_trans:

yn_out_dns_trans
^^^^^^^^^^^^^^^^

[currently not available] character default='n'

   | Available for ``theory='tddft_pulse'``.
   | Switch to calculate transition density at specified frequency omega (specified by ``out_dns_trans_energy``), drho(r,omega)=FT(rho(r,t)-rho_gs(r))/T.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_dns_trans_energy:

out_dns_trans_energy
^^^^^^^^^^^^^^^^^^^^

[currently not available] real(8), default=1.55d0 [eV]

   | Available for ``theory='tddft_pulse'``.
   | A frequency to output drho(r,omega)=FT(rho(r,t)-rho_gs(r))/T.

.. _yn_out_elf:

yn_out_elf
^^^^^^^^^^

character, default='n'

   | Available for ``theory='dft'``.
   | Switch to output the electron localization function.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_out_elf_rt:

yn_out_elf_rt
^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='dft_md', 'tddft_pulse'``.
   | Switch to output the electron localization function during the time propagation every ``out_elf_rt_step`` time steps.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_elf_rt_step:

out_elf_rt_step
^^^^^^^^^^^^^^^

integer, default=50

   | Available for ``theory='dft_md', 'tddft_pulse'``.
   | Electron localization function during the time propagation is outputted every ``out_elf_rt_step`` time steps.

.. _yn_out_estatic_rt:

yn_out_estatic_rt
^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='tddft_pulse'``.
   | Switch to print the static electric field during the time propagation every ``out_estatic_rt_step`` time steps.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_estatic_rt_step:

out_estatic_rt_step
^^^^^^^^^^^^^^^^^^^

integer, default=50

   | Available for ``theory='tddft_pulse'``.
   | The static electric field during the time propagation is outputed every ``out_estatic_rt_step`` time steps.

.. _yn_out_rvf_rt:

yn_out_rvf_rt
^^^^^^^^^^^^^

character, default='n'

   | Available for TDDFT based options and 'dft_md' option of ``theory``.
   | Switch to print the coordinates[A], velocities[au], forces[au] of atoms during time-propagation in ``SYSname``\_trj.xyz every ``out_rvf_rt_step`` time steps.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable
   | If ``yn_md='y'``, this option is automatically turned on.

.. _out_rvf_rt_step:

out_rvf_rt_step
^^^^^^^^^^^^^^^

integer, default=10

   | Available for TDDFT based options and 'dft_md' option of ``theory``.
   | The coordinates[A], velocities[au], forces[au] of atoms during time-propagation are outputed in ``SYSname``\_trj.xyz every ``out_rvf_rt_step`` time steps.

.. _yn_out_tm:

yn_out_tm
^^^^^^^^^

[Trial] character, default='n'

   | Available for ``yn_periodic='y'`` with ``theory='dft'``.
   | Switch to calculate and print the transition matrix elements between occupied and virtual orbitals to ``SYSname``\_tm.data after the ground state calculation.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_out_gs_sgm_eps:

yn_out_gs_sgm_eps
^^^^^^^^^^^^^^^^^

[Trial] character, default='n'

   | Available for ``theory='dft'``.
   | Switch to calculate and print conductivity (sigma) and dielectric function (epsilon) based on transition moment after convergence of the ground state calculation. These are printed in the output files, ``SYSname``\_sigma.data and ``SYSname``\_epsilon.data
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_gs_sgm_eps_mu_nu

out_gs_sgm_eps_mu_nu
^^^^^^^^^^^^^^^^^^^^

integer, default=3,3

   | Available for ``yn_out_gs_sgm_eps='y'`` with ``theory='dft'``.
   | Index of conductibity and dielectric tensol element calculated in this option. Default of (3,3) means zz element.

.. _out_gs_sgm_eps_width

out_gs_sgm_eps_width
^^^^^^^^^^^^^^^^^^^^

real(8), default=0.015d0 [eV]

   | Available for ``yn_out_gs_sgm_eps='y'`` with ``theory='dft'``.
   | Smearing width used in conductivity and dielectric function

.. _out_ms_step:

out_ms_step
^^^^^^^^^^^

integer, default=100

   | Available for ``theory='multi_scale_maxwell_tddft'``.
   | Some quantities are printed every ``out_ms_step`` time step in the Maxwell-TDDFT multiscale calculations.

.. _format_voxel_data:

format_voxel_data
^^^^^^^^^^^^^^^^^

character, default='cube'

   | Available for ``yn_out_psi='y'``, ``yn_out_dns(_rt)='y'``,  ``yn_out_dns_ac_je='y'``,  ``yn_out_elf(_rt)='y'``,  ``yn_out_estatic_rt='y'``.
   | Option of the file format for three-dimensional volumetric data.
   |   ``'avs'`` /  AVS format
   |   ``'cube'`` / cube format
   |   ``'vtk'`` / vtk format

.. _nsplit_voxel_data:

nsplit_voxel_data
^^^^^^^^^^^^^^^^^

integer, default=1

   | Available for ``format_voxel_data='avs'``.
   | Number of separated files for three dimensional data.

.. _yn_lr_w0_correction:

yn_lr_w0_correction
^^^^^^^^^^^^^^^^^^^

[Trial] character, default='n'

   | Available for ``yn_periodic='y'`` and ``trans_longi='tr'`` with ``theory='tddft_response'``.
   | Apply correction around zero frequency of dielectric function to suppress numerical error.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_out_current_decomposed:

yn_out_current_decomposed
^^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for TDDFT based options of ``theory``.
   | Switch to output docomposed elements of the electron current density.
   | The sum of the docomposed elements is equal to the current density in SYSname_rt.data.
   | Output file: SYSname_current_decomposed.data
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _out_current_decomposed_step:

out_current_decomposed_step
^^^^^^^^^^^^^^^^^^^^^^^^^^^

integer, default=100

   | Available when ``yn_out_current_decomposed='y'``.
   | The decomposed current data is outputted every ``out_current_decomposed_step`` step.
   
.. _out_rt_spin_step:

out_rt_spin_step
^^^^^^^^^^^^^^^^

integer, default=100

   | Available for TDDFT based methods with ``spin='noncollinear'``.
   | The spin magnetization and spin current density are outputted every ``out_rt_spin_step`` time steps in the output file SYSname_rt_spin.data.
   | For the definition of the spin current, see [N. Tancogne-Dejean et al, npj Computational Materials 8, 145 (2022).].
   
.. _yn_out_mag_decomposed_rt:

yn_out_mag_decomposed_rt
^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for TDDFT based methods with ``spin='noncollinear'``.
   | Switch to output docomposed elements of the time-dependent spin magnetization at every ``out_rt_spin_step`` time steps.
   | Output file: SYSname_mag_decomposed_rt.data
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable
   
.. _yn_out_spin_current_decomposed:

yn_out_spin_current_decomposed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for TDDFT based methods with ``spin='noncollinear'``.
   | Switch to output docomposed elements of the spin current density at every ``out_rt_spin_step`` time steps.
   | Output file: SYSname_spin_current_decomposed.data
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable
   
.. _yn_out_spin_current_micro:

yn_out_spin_current_micro
^^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for TDDFT based methods with ``spin='noncollinear'``.
   | Switch to output voxel data files of the microscopic spin-current density at every ``out_rt_spin_step`` time steps.
   | Output files: spin_curr_micro_[xyz]_[xyz]_000001.<format_voxel_data> etc.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _yn_out_perflog:

yn_out_perflog
^^^^^^^^^^^^^^

character, default='y'

   | Available for all ``theory``
   | Switch to print the performance log of routines and modules.
   | Options:
   |   ``'y'`` / enable
   |   ``'n'`` / disable

.. _format_perflog:

format_perflog
^^^^^^^^^^^^^^

character, default='stdout'

   | Available for ``yn_out_perflog = 'y'``
   | The output format of performance log.
   | Options:
   |   ``'stdout'`` / standard output unit
   |   ``'text'`` / save as a text file
   |   ``'csv'`` / save as a csv format file

.. _&poisson:

&poisson
--------

.. _method_poisson:

method_poisson
^^^^^^^^^^^^^^^^

character, Default='cg'

   | Available for ``yn_periodic='n'`` in DFT and TDDFT based options of ``theory``.
   | This papameter specify how to solve the Poisson equation.
   | Options:
   |  ``cg``/ Conjugate-Gradient(CG) method
   |  ``ft``/ Fourier transformation method
   |  ``dirichlet``/ Dirichlet boundary condition method

.. _layout_multipole:

layout_multipole
^^^^^^^^^^^^^^^^

character, Default=3

   | Available for ``yn_periodic='n'`` in DFT and TDDFT based options of ``theory``.
   | This papameter specify how to achieve multipole expansioin in the Hartree potential calculation.
   | Options:
   |  ``1``/ A single pole at the center.
   |  ``2``/ Multipoles are set at each center of atoms.
   |  ``3``/ Multipoles are set at the center of mass of electrons in prepared cuboids in each process.

.. _num_multipole_xyz(3):

num_multipole_xyz(3)
^^^^^^^^^^^^^^^^^^^^

integer, default=0

   | Available for ``yn_periodic='n'`` in DFT and TDDFT based options of ``theory``.
   | Number of multipoles. When default is set, the number of multipoles is calculated automatically.

.. _lmax_multipole:

lmax_multipole
^^^^^^^^^^^^^^

[Trial] integer, default=4

   | Available for ``yn_periodic='n'`` in DFT and TDDFT based options of ``theory``.
   | A maximum order of the multipole expansion to prepare boundary condition of Poisson equation.

.. _threshold_cg:

threshold_cg
^^^^^^^^^^^^

real(8), default=1d-15 [a.u.]

   | Available for ``yn_periodic='n'`` in DFT and TDDFT based options of ``theory``.
   | A threshold for the convergence of the Hartree-cg calculation. A quantity examined is given by ||tVh(i)-tVh(i-1)||^2/(number of grids).

.. _&ewald:

&ewald
------

.. _newald:

newald
^^^^^^

integer, default=4

   | Available for ``yn_periodic='y'`` in DFT/TDDFT based options of ``theory``.
   | Parameter of the Ewald method for the ion-ion Coulombic interaction. Short-range part of the Ewald sum is calculated within ``newald``-th nearlist neighbor cells.

.. _aewald:

aewald
^^^^^^

real(8), default=0.5d0 [a.u.]

   | Available for ``yn_periodic='y'`` in DFT/TDDFT based options of ``theory``.
   | Square of range separation parameter for Ewald method (This parameter is given only in atomic unit).

.. _cutoff_r:

cutoff_r
^^^^^^^^

real(8), default=-1d0

   | Available for ``yn_periodic='y'`` in DFT/TDDFT based options of ``theory``.
   | Cut-off length in real-space. The length is automatically determined if ``cutoff_r`` < 0.

.. _cutoff_r_buff:

cutoff_r_buff
^^^^^^^^^^^^^

real(8), default=2d0 [a.u.]

   | Available for ``yn_periodic='y'`` in ``yn_md='y'`` or in ``theory='dft_md'``.
   | Buffer length in radius for book-keeping for real-space interaction.

.. _cutoff_g:

cutoff_g
^^^^^^^^

real(8), Default=-1d0

   | Available for ``yn_periodic='y'`` in DFT/TDDFT based options of ``thddeory``.
   | Cut-off in G-space in the Ewald method. No cut-off in default.

.. _&opt[Trial]:

&opt[Trial]
-----------

.. _nopt:

nopt
^^^^

integer, default=100

   | Available for ``yn_opt='y'`` in ``theory='dft'``.
   | The maximum step number of geometry optimization.

.. _convrg_opt_fmax:

convrg_opt_fmax
^^^^^^^^^^^^^^^

real(8), default=1d-3 (a.u.)

   | Available for ``yn_opt='y'`` in ``theory='dft'``.
   | Convergence threshold of geometry optimization is specified for the maximum force acting on atoms.

.. _max_step_len_adjust:

max_step_len_adjust
^^^^^^^^^^^^^^^^^^^

real(8), default=-1d0

   | Available for ``yn_opt='y'`` in ``theory='dft'``.
   | Set maximum optimization step length (if positive number is given)

.. _&md[Trial]:

&md[Trial]
-----------

.. _ensemble:

ensemble
^^^^^^^^

character, default='NVE'

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Ensemble in MD option:
   | Options:
   |  ``NVE``/ NVE ensemble (constant energy and volume system)
   |  ``NVT``/ NVT ensemble (constant temperature and volume system)

.. _thermostat:

thermostat
^^^^^^^^^^

character, default='nose-hoover'

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Thermostat in "NVT" option:
   | Options:
   |  ``nose-hoover``/ Nose-Hoover thermostat

.. _step_velocity_scaling:

step_velocity_scaling
^^^^^^^^^^^^^^^^^^^^^

integer, default=-1

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Time step interval for velocity-scaling. Velocity-scaling is applied if this is set to positive.

.. _step_update_ps:

step_update_ps
^^^^^^^^^^^^^^

integer, default=10

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Time step interval for updating pseudopotential (Larger number reduces computational time but increases inaccuracy).

.. _temperature0_ion_k:

temperature0_ion_k
^^^^^^^^^^^^^^^^^^

real(8), Default=298.15d0 [K]

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Setting ionic temperature in unit of [K] for NVT ensemble, velocity scaling and generating initial velocities.

.. _yn_set_ini_velocity:

yn_set_ini_velocity
^^^^^^^^^^^^^^^^^^^

character, Default='n'

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Switch to generate initial velocities.
   | Options:
   |  ``y``/ Generate initial velocity with Maxwell-Bortzman distribution
   |  ``n``/ disable

.. _file_ini_velocity:

file_ini_velocity
^^^^^^^^^^^^^^^^^

[Trial] character, default='none'

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | File name for reading initial velocities. This is read if the file name is given, then, the priority is higher than use of ``set_ini_velocity`` and restart data of velocities. The format is simply vx(iatom) vy(iatom) vz(iatom) in each line. The order of atoms must be the same as the given coordinates in the main input file. In case of using nose-hoover thermostat, a thermostat variable should be put at the last line (all atomic unit).

.. _thermostat_tau:

thermostat_tau
^^^^^^^^^^^^^^

real(8), default=1d0 [fs]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Parameter in Nose-Hoover method: controlling time constant for temperature.

..
   #XXX removed?#
   - **seed_ini_velocity** (integer, Default=123)
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Random seed (integer number) to generate initial velocity if ``set_ini_velocity`` is set to y.
   Default is ``123``.
..

.. _yn_stop_system_mom:

yn_stop_system_mom
^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   | Center of mass is fixed every time step.
   | Options:
   |  ``y``/ enable
   |  ``n``/ disable

.. _&jellium:

&jellium
--------

.. _yn_jm:

yn_jm
^^^^^

character, default='n'

   | Available for the DFT/TDDFT based options of ``theory``.
   | Switch to use jellium model.
   | Options:
   |  ``y``/ enable
   |  ``n``/ disable
   | When ``yn_jm='y'``, ``&functional/xc`` must be ``'pz'``.

.. _yn_charge_neutral_jm:

yn_charge_neutral_jm
^^^^^^^^^^^^^^^^^^^^

character, default='y'

   | Available for ``yn_jm='y'`` in the DFT/TDDFT based options of ``theory``.
   | Option to enforce exact charge neutrality :
   | Options:
   |  ``y``/ enable. ``rs_bohr_jm`` is modified to fulfill exact charge neutrality.
   |  ``n``/ disable. ``rs_bohr_jm`` is not modified, and there may appears small charge-neutrality error.

.. _yn_output_dns_jm:

yn_output_dns_jm
^^^^^^^^^^^^^^^^

character, default='y'

   | Available for ``yn_jm='y'`` in the DFT/TDDFT based options of ``theory``.
   | Switch to output positive background charge density.
   | Options:
   |  ``y``/ enable
   |  ``n``/ disable

.. _shape_file_jm:

shape_file_jm
^^^^^^^^^^^^^

character, default='none'

   | Available for ``yn_jm='y'`` in the DFT/TDDFT based options of ``theory``.
   | Name of input shape file that contains positive background charge density to be used in the jellium model calculations. The shape file can be generated by using ``FDTD_make_shape`` in SALMON utilities (https://salmon-tddft.jp/utilities.html). When ``shape_file_jm='none'``, the shape of the positive background charge density is specified by ``sphere_nion_jm`` and ``sphere_loc_jm`` which generate spherical shapes.

.. _num_jm:

num_jm
^^^^^^

integer, Default=0

   | Available for ``yn_jm='y'`` in the DFT/TDDFT based options of ``theory``.
   | When ``shape_file_jm`` is not 'none', ``num_jm`` specifies number of media used in the jellium model. When ``shape_file_jm='none'``, ``num_jm`` specifies number of spherical shapes.

.. _rs_bohr_jm(:):

rs_bohr_jm(:)
^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``yn_jm='y'`` in the DFT/TDDFT based options of ``theory``.
   | When ``shape_file_jm`` is not 'none', ``rs_bohr_jm(n)`` spesifies the Wigner-Seitz radius of n-th media. When ``shape_file_jm='none'``, ``rs_bohr_jm(n)`` spesifies the Wigner-Seitz radius of n-th sphere.

.. _sphere_nelec_jm(:):

sphere_nion_jm(:)
^^^^^^^^^^^^^^^^^^

integer, default=0

   | Available for ``yn_jm='y'`` and ``shape_file_jm='none'`` in the DFT/TDDFT based options of ``theory``. ``sphere_nion_jm(n)`` spesifies ion number for n-th sphere. At present, only neutral systems can be treated.

.. _sphere_loc_jm(:,3):

sphere_loc_jm(:,3)
^^^^^^^^^^^^^^^^^^

real(8), default=0d0

   | Available for ``yn_jm='y'`` and ``shape_file_jm='none'`` in the DFT/TDDFT based options of ``theory``. ``sphere_loc_jm(n,1:3)=x,y,z`` spesifies location of center of mass for n-th sphere. Note that the coordinate system ranges from ``-al/2`` to ``al/2`` for ``&system/yn_periodic='n'`` while ranges from ``0`` to ``al`` for ``&system/yn_periodic='y'``.

.. _&code:

&code
-----

.. _yn_want_stencil_hand_vectorization:

yn_want_stencil_hand_vectorization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

character, default='y'

   | Switch to use hand-vectorized optimization code of stencil in the hamiltonian calculation.
   | SALMON checks if the calculation can use the hand-vectorized code. If it fails, SALMON will use a typical implementation.

.. _yn_want_communication_overlapping:

yn_want_communication_overlapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

character, default='n'

   | Available for ``theory='tddft*' or '*maxwell_tddft'``
   | Switch to use computation/communication overlap algorithm to improve the performance of stencil in the hamiltonian calculation.
   | SALMON checks if the calculation can use the overlap algorithm. If it fails, SALMON will uses a non-overlap algorithm.

.. _stencil_openmp_mode:

stencil_openmp_mode
^^^^^^^^^^^^^^^^^^^

character, default='auto'

   | This option selects an OpenMP parallelization mode of stencil in the hamiltonian calculation.
   | Options:
   |   ``auto``    / SALMON decides the parallelization target automatically.
   |   ``orbital`` / OpenMP parallelization is applied to orbital (and k-point) loop.
   |   ``rgrid``   / OpenMP parallelization is applied to real-space grid loop.

.. _current_openmp_mode:

current_openmp_mode
^^^^^^^^^^^^^^^^^^^

character, default='auto'

   | This selects an OpenMP parallelization mode of the current calculation.
   | Options:
   |   ``auto``    / SALMON decides the parallelization target automatically.
   |   ``orbital`` / OpenMP parallelization is applied to orbital (and k-point) loop.
   |   ``rgrid``   / OpenMP parallelization is applied to real-space grid loop.

.. _force_openmp_mode:

force_openmp_mode
^^^^^^^^^^^^^^^^^

character, default='auto'

   | This selects an OpenMP parallelization mode of the force calculation.
   | Options:
   |   ``auto``    / SALMON decides the parallelization target automatically.
   |   ``orbital`` / OpenMP parallelization is applied to orbital (and k-point) loop.
   |   ``rgrid``   / OpenMP parallelization is applied to real-space grid loop.

..
  #### Following keywords are commented out as these are originated from GCEED and supposed to be removed ####
  **Following variables are moved from the isolated part. Some of them may be added to common input, be combined to it, and be removed.**

  &group_fundamental[Trial]
  -------------------------

  - **iwrite_projection** (integer, Default=0)[Trial]
   | Available for ``theory='XXX'``.
   A variable for projection.

  - **itwproj** (integer, Default=-1)[Trial]
   | Available for ``theory='XXX'``.
   The projection is calculated every ``itwproj`` time steps.

  - **iwrite_projnum** (integer, Default=0)[Trial]
   | Available for ``theory='XXX'``.
   There is a malfunction in this variable.

  &group_others[Trial]
  ---------------------

  - **num_projection** (Interger, Default=1)[Trial]
   | Available for ``theory='XXX'``.
   |  of orbitals for projections.

  - **iwrite_projection_ob(200)** (Interger, Default=1, 2, 3, ..., 200)[Trial]
   | Available for ``theory='XXX'``.
   Orbital number to be written as projections.

  - **iwrite_projection_k(200)** (Interger, Default=1)[Trial]
   | Available for ``theory='XXX'``.
   This variable will be removed.
