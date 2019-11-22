.. _List of all input keywords:

List of all input keywords
==========================

-  `&calculation`_
-  `&control`_
-  `&units`_
-  `&parallel`_
-  `&system`_
-  `&atomic_red_coor`_
-  `&atomic_coor`_
-  `&pseudo`_
-  `&functional`_
-  `&rgrid`_
-  `&kgrid`_
-  `&tgrid`_
-  `&propagation`_
-  `&scf`_
-  `&emfield`_
-  `&multiscale`_
-  `&maxwell`_
-  `&analysis`_
-  `&poisson`_
-  `&ewald`_
-  `&opt(Trial)`_
-  `&md(Trial)`_
-  `&group_fundamental(Trial)`_
-  `&group_parallel(Trial)`_  
-  `&group_hartree(Trial)`_ 
-  `&group_file(Trial)`_
-  `&group_others(Trial)`_


&calculation
------------

- **theory** (character, 0d/3d)
   | Choice of Calculation theories.
   |  Options
   |    ``DFT``  / ground state calculation based on DFT
   |    ``DFT_MD``  / adiabatic ab initio MD simulations based on DFT
   |    ``TDDFT_response``  / simulations based on TDDFT for response
   |    ``TDDFT_pulse``  / simulations based on TDDFT using pulsed light
   |    ``Single_scale_Maxwell_TDDFT``  / coupled Maxwell and TDDFT single-scale simulation
   |    ``Multi_scale_Maxwell_TDDFT``  / coupled Maxwell and TDDFT multi-scale simulation
   |    ``Maxwell``  / electromagnetic simulations based on the Maxwell's equations (o

- (Trial) **yn_md** (character, 0d/3d)
   Enable(``'y'``)/disable(``'n'``). 
   Molecular dynamics option. This is available for ``theory='DFT'`` (Adiabatic ground-state MD) and ``theory='TDDFT_pulse'`` (Ehrenfest MD)
   Default is ``'n'``.

- (Trial) **yn_opt** (character, 0d/3d)
   Enable(``'y'``)/disable(``'n'``). 
   Geometry optimization option. Available for ``theory='DFT'``
   Default is ``'n'``.


&control
--------

- **sysname** (character, 0d/3d)
   Name of calculation. This is used for a prefix of output files.
   Default is ``default``.

- **base_directory** (character, 0d/3d)
   Name of a default directory, where the basic results will be written down.
   Default is the current directoy, ``./``.

- **output_buffer_interval** (integer, 0d/3d)
   xxx.

- (Trial) **yn_restart** (character, 3d)
   Flag for restart, ``'new'`` or ``'restart'``.
   ``'new'`` is default.

- (Trial) **directory_read_data** (character, 3d)
   xxx.

- (Trial) **yn_self_checkpoint** (character, 3d)
   xxx.

- (Trial) **checkpoint_interval** (integer, 3d)
   Frequency of backup during the time-propagation. 
   If ``0`` is set, the backup is not performed.
   Default is ``0``.

- (Trial) **time_shutdown** (real(8), 3d)
   Timer for automatic shutdown. The unit is always second.
   If negative time is chosen, the automatic shutdown is not performed.
   Default is ``-1`` sec.


&units
------

- **unit_system** (character, 0d/3d)
   Unit for input variables. 
   If ``'au'`` or ``'a.u.'``, atomic unit system is used. 
   If ``'A_eV_fs'``, Angstrom-eV-fs unit system is used. 


&parallel
---------

- (Trial) **yn_domain_parallel** (character, 3d)
   If specified ``yn_domain_parallel='y'`` and ``&system/iperiodic=3``, program codes for domain parallel version run in periodic system calculations.

- **nproc_k/nproc_ob/nproc_domain_orbital(3)/nproc_domain_general(3)** (integer, 0d)
   Followings are explanation of each variable.

  - ``nproc_k``: Number of MPI parallelization for orbitals that related to the wavefunction calculation.
  - ``nproc_ob``: Number of MPI parallelization for orbitals that related to the wavefunction calculation.
  - ``nproc_domain_orbital(3)'``: Number of MPI parallelization for each direction in real-space that related to the wavefunction calculation. 
  - ``nproc_domain_general(3)'``: Number of MPI parallelization for each direction in real-space that related to the electron density calculation. 

    Defaults are ``0`` for ``nproc_k``/``nproc_ob`` and ``(0/0/0)`` for ``nproc_domain_orbital``/``nproc_domain_s``. If users use the defaults, automatic proccess assignment is done. Users can also specify ``nproc_k``, ``nproc_ob``, ``nproc_domain``, and ``nproc_domain_general`` manually. In that case, ``nproc_k`` must be set to ``1`` for isolated system calculations. In addition, followings must be satisfied.

  - ``nproc_k`` \* ``nproc_ob`` \* ``nproc_domain_orbital(1)`` \* ``nproc_domain_orbital(2)`` \* ``nproc_domain_orbital(3)`` \= total number of processors
  - ``nproc_domain_general(1)`` \* ``nproc_domain_general(2)`` \* ``nproc_domain_general(3)`` \= total number of processors
  - ``nproc_domain_general(1)`` is a multiple of ``nproc_domain_orbital(1)``
  - ``nproc_domain_general(2)`` is a multiple of ``nproc_domain_orbital(2)``
  - ``nproc_domain_general(3)`` is a multiple of ``nproc_domain_orbital(3)``

- **num_datafiles_in/num_datafiles_out** (integer, 0d)
   Number of input/output files for wavefunction.
   Defaults are ``1``. If ``num_datafiles_in``/``num_datafiles_out`` are 1, wave functions are read from/ written in a regular intermediate file. If ``num_datafiles_in``/``num_datafiles_out`` are larger than or equal to 2, the wave functions are read from/ written in separated intermediate data files, and number of files are equal to ``num_datafiles_in``/``num_datafiles_out``. These variables must be equal to nth power of 2. (n: 0 or positive integer)

- **yn_ffte** (character, 0d)
   Method of Fourier transformation.  ``'ft'``,  ``'FT'``, ``'ffte'`` or ``'FFTE'`` can be chosen.
   Default is ``'ft'``.
   This variable is effective only when ``yn_domain_parallel='y'`` and ``&system/iperiodic=3``.

- **process_allocation** (character, 0d)
   xxx.


&system 
-------

- **yn_periodic** (character, 0d/3d)
   Dimension for periodic boundary condition.
   ``0`` is for isolated systems, and 
   ``3`` is for solids.
   Default is ``0``.

- **ispin** (integer, 0d)
   Variable for classification of closed shell systems and open shell systems.
   ``0`` is for closed shell systems, and
   ``1`` is for open shell systems.
   Default is ``0``

- **al(3)** (real(8), 0d/3d)
   Lattice constants. Unit of the length can be chosen by ``&units/unit_system``.

- **al_vec1(3)/al_vec2(3)/al_vec3(3)** (real(8), 3d)
   xxx.

- **isym** (integer, 3d)
   Number of symmetries that can be used for reduction of k-points.
   Default is ``0``.

- **crystal_structure** (character, 3d)
   Name of symmetry that can be used for the reduction of # of k-points.
   Default is ``'none'``.

- **nstate** (integer, 0d/3d)
   Number of states/bands.

- **nstate_spin(2)** (integer, 0d)
   Number of states/bands can be specified independently by ``nstate_spin(1)/nstate_spin(2)``.
   This option is incompatible with ``nstate``

- **nelec** (integer, 0d/3d)
   Number of valence electrons.

- **nelec_spin(2)** (integer, 0d)
   Number of up/down-spin electrons can be specified independently by ``nelec_spin(1)/nelec_spin(2)``.
   This option is incompatible with ``nelec``

- **temperature** (real(8), 3d)
   Temperature of electrons. When you calculate a system of zero band-gap energy like metals, zero or positive number of the temperature should be given.
   Unit of the energy can be chosen ``&units/unit_system``. 
   Default is ``-1.0`` (this is for system which has a band gap energy).

- (Trial) **temperature_k** (real(8), 0d)
   Temperature of electrons [K]. Default is ``-1.d0``.

- **nelem** (integer, 0d/3d)
   Number of elements that will be used in calculations.

- **natom** (integer, 0d/3d)
   Number of atoms in a calculation cell.


- (Trial) **file_atom_red_coor** (character, 3d)
   File name of atomic positions. In this file, 
   the atomic coordinates can be written in reduced coordinates.
   This option is incompatible with 
   ``&system/file_atom_coor``,
   ``&atomic_coor``, and 
   ``&atomic_red_coor``.

- (Trial) **file_atom_coor** (character, 0d)
   File name of atomic positions. In this file, 
   the atomic coordinates can be written in Cartesian cooridnates.
   The unit of the length can be chosen by 
   ``&units/unit_system``.
   This option is incompatible with 
   ``&system/file_atom_red_coor``,
   ``&atomic_coor``, and 
   ``&atomic_red_coor``.


&atomic_red_coor
----------------

In ``&atomic_red_coor``, positions of atoms can be written in reduced coordinates
as follows:

|  'Si'	 0.00  0.00  0.00  1
|  'Si'	 0.25  0.25  0.25  1
|  ...

Here, the information of atoms is ordered in row. For example, the first row gives
the information of the first atom. The number of rows must be equal to 
``&system/natom``.
The first coloum can be any caracters and does not affect calculations.
The second, third and fourth columns are reduced coordinates for
the first, second and third directions, respectively. 
The fifth column is a serial number of the atom spieces, which is used in 
``&pseudo``.
This option is incompatible with 
``&system/file_atom_red_coor``,
``&system/file_atom_coor``, and
``&atomic_coor``.


&atomic_coor
------------

In &atomic_coor, positions of atoms can be written in Cartesian coordinates.
The structure is same as &atomic_red_coor.
The unit of the length can be chosen by 
``&units/unit_length``.
This option is incompatible with 
``&system/file_atom_red_coor``,
``&system/file_atom_coor``, and
``&atomic_red_coor``.


&pseudo
-------

Input for psudopotentials. Size of array (:) is equal to ``&system/nelem``.

- **file_pseudo(:)** (character, 0d/3d)
   Name of pseudopotential files.

- **lmax_ps(:)** (integer, 0d/3d)
   Maximum angular momentum of pseudopotential projectors.

- **lloc_ps(:)** (integer, 0d/3d)
   Angular momentum of pseudopotential that will be treated as local.

- **izatom(:)** (integer, 0d/3d)
   Atomic number.

- (Trial) **yn_psmask(:)** (character, 0d/3d)
   Enable(``'y'``)/disable(``'n'``) 
   Fourier filtering for pseudopotentials. 
   Default is ``'n'``.

- (Trial) **alpha_mask(:)** (real(8), 0d/3d)
   Parameter for the Fourier filtering for pseudopotential.
   Default is ``'0.8'``.

- (Trial) **gamma_mask(:)** (real(8), 0d/3d)
   Parameter for the Fourier filtering for pseudopotential.
   Default is ``'1.8'``.

- (Trial) **eta_mask(:)** ``Real(8)``); 0d/3d)
   Parameter for the Fourier filtering for pseudopotential.
   Default is ``'15.0'``.


&functional
-----------

- **xc** (character, 0d/3d)
   Exchange-correlation functionals.
   At present version, the functional 'PZ', 'PZM' and 'TBmBJ' is available for both 0d/3d calculations, and the functionals 'TPSS' and 'VS98' are available for 3d calculations.

  - ``'PZ'``: Perdew-Zunger LDA :Phys. Rev. B 23, 5048 (1981).
  - ``'PZM'``: Perdew-Zunger LDA with modification to improve sooth connection between high density form and low density one. :J. P. Perdew and Alex Zunger, Phys. Rev. B 23, 5048 (1981).
  - ``'TBmBJ'``: Tran-Blaha meta-GGA exchange with Perdew-Wang correlation. :Fabien Tran and Peter Blaha, Phys. Rev. Lett. 102, 226401 (2008). John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992).
  - ``'TPSS'``: Tao, Perdew, Staroverov and Scuseria meta-GGA exchange correlation. :J. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria, Phys. Rev. Lett. 91, 146401 (2003).
  - ``'VS98'``:  van Voorhis and Scuseria exchange with Perdew-Wang correlation: T. Van Voorhis and G. E. Scuseria, J. Chem. Phys. 109, 400 (1998).

- **cname, xname** (character, 0d/3d)
   xxx.

- **alibxc, alibx, alibc** (character, 0d/3d)
   By specifying ``alibxc``, the functionals prepared in libxc package are available. 
   They can be set indivisually by specifying ``alibx`` and ``alibc``.
   To use libxc libraries, ``--with-libxc`` option must be added in excecuting configure. 
   The available option of the exchange-correlation functionals are listed in the LibXC website. 
   [See http://www.tddft.org/programs/libxc/functionals/]
   
- **cval** (real(8), 3d)
   Mixing parameter in Tran-Blaha meta-GGA exchange potential. If ``cval`` is set to a minus value, the mixing-parameter computed
   by the formula in the original paper [Phys. Rev. Lett. 102, 226401 (2008)].
   Default is estimated from :math:`\left\langle |\nabla \rho(\mathbf{r};t)| / \rho(\mathbf{r};t) \right\rangle`.


&rgrid
------

- **dl(3)** (real(8), 0d/3d)
   Spacing of real-space grids. Unit of length can be chosen by
   ``&units/unit_system``.
   This valiable cannot be set with 
   ``&rgrid/num_rgrid`` simultaneously.
   If ``&system/iperiodic`` is set to ``3``,
   the actual grid spacing is automatically refined in calculations
   so that the size of the simulation box
   ``&system/al(3)`` becomes divisible by the spacing.

- **num_rgrid(3)** (integer, 3d)
   Number of real-space grids.
   This valiable cannot be set with 
   ``&rgrid/dl`` simultaneously.


&kgrid
------

- **num_kgrid(3)** (integer, 3d)
   Number of k-points (grid points of k-vector) discretizing
   the Brillouin zone.
   Each component must be even.

- **file_kw** (character, 3d)
   Name of a file for flexible k-point sampling.
   This file will be read if ``num_kgrid`` is smaller than 1.


&tgrid
------

- **nt** (integer, 0d/3d)
   Number of total time steps for real-time propagation.

- **dt** (real(8), 0d/3d)
   Time step. Unit of time can be chosen by ``&units/unit_system``.


&propagation
------------

- **propagator** (character, 3d)
   Choice of Propagator.
   ``middlepoint`` is an propagator
   with the Hamiltoinan at midpoint of two-times.
   ``etrs`` is enforced time-reversal symmetry propagator.
   [M.A.L. Marques, A. Castro, G.F. Bertsch, and A. Rubio, Comput. Phys. Commun., 151 60 (2003)].
   Default is ``middlepoint``.

- (Trial) **n_hamil** (integer, 0d)
   Order of Taylor expansion of a propagation operator.
   Default is ``4``.

- (Trial) **yn_fix_func** ``character(1)``; 3d)
   Option not to update functional (or Hamiltonian) in RT calculation, i.e., keep ground state Hamiltonian during time-evolution.
   Default is ``'n'``.

&scf
----

- **nscf** (integer, 0d/3d)
   Number of maximum scf cycle.

- **ncg** (integer, 0d/3d)
   Number of interation of Conjugate-Gradient method for each scf-cycle.
   Default is ``5``.

- **method_mixing** (character, 0d) 
   Methods for density/potential mixing for scf cycle. ``simple`` and ``broyden`` can be chosen.
   Default is ``broyden``.

- **mixrate** (real(8), 0d)
   Mixing ratio for simple mixing. Default is ``0.5``.

- **nmemory_mb** (integer, 0d/3d)
   Number of stored densities at previous scf-cycles for 
   the modified-Broyden method. Default is ``8``. 
   If ``&system/iperiodic`` is ``0``, ``nmemory_mb`` must be less than 21.

- **alpha_mb** (real(8), 0d/3d)
   Parameter of the modified-Broyden method.
   Default is ``0.75``.

- **yn_subspace_diagonalization** (character, 0d)
   Enable(``'y'``)/disable(``'n'``) 
   subspace diagonalization during scf cycle.

- **convergence** (character, 0d/3d)
   Choice of quantity that is used for convergence check in a scf calculation. 
   Default is ``'rho_dne'``. 

  - ``'rho_dne'``: Convergence is checked by sum_ix|rho(ix,iter)-rho(ix,iter-1)|dx/N, where iter is an iteration number of the scf calculation and N is ``&system/nelec``, the number of the valence electrons.

   For isolated systems, the followings can also be chosen.

  - ``'norm_rho'``: Convergence is checked by the square of the norm of difference of density, ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\=sum_ix|rho(ix,iter)-rho(ix,iter-1)|\ :sup:`2`\. 
  - ``'norm_rho_dng'``: Convergence is checked by ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\/(number of grids). "dng" means "devided by number of grids".
  - ``'norm_pot'``: Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\, where Vlocal is Vh + Vxc + Vps_local.
  - ``'pot_dng'``: Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\/(number of grids).

- **threshold** (real(8), 0d/3d)
   Threshold for convergence check that is used when ``'rho_dne'`` is specified.
   Default is ``1d-17``. 
   XXX(threshold_norm_rho (real(8), 0d))XXX
   Threshold for convergence check that is used when either ``'norm_rho'`` or ``'norm_rho_dng'`` is specified. ``threshold_norm_rho`` must be set when either ``'norm_rho'`` or ``'norm_rho_dng'`` is specified.
   Default is ``-1d0`` a.u. (1 a.u.= 45.54 A\ :sup:`-6`\)
   XXX(threshold_norm_pot (real(8), 0d))XXX
   Threshold for convergence check that is used when either ``'norm_pot'`` or ``'norm_pot_dng'`` is specified. ``threshold_norm_pot`` must be set when either ``'norm_pot'`` or ``'norm_pot_dng'`` is specified.
   Default is ``-1d0`` a.u. (1 a.u.= 33.72x10\ :sup:`4`\ A\ :sup:`-6`\eV\ :sup:`2`\)

- **omp_loop** (character, 3d)
   XXX only ARTED XXX
   Loop for OpenMP parallelization in the ground state SCF if periodic boundary system is used. 

  - ``k``: parallelization for k-point loop (Default).
  - ``b``: parallelization mainly for band orbital loop (sometimes space grid loop too). This works efficiently if the number of k-point treated in each node is small (e.x. the case of single k-point for each node)


- (Trial) **skip_gsortho** (character, 3d)
   XXX only ARTED XXX
   Flag to skip Gram-Schmidt orthogonalization in CG loop if periodic boundary system is used. If this is skipped the more iteration number is necessary to get convergence but each iteration step gets faster. If ``omp_loop=b``, this flag is always applied.
   Default is ``n``



&emfield
--------

- **trans_longi** (character, 3d)
   Geometry of solid-state calculations.
   Transverse ``'tr'`` and longitudinal ``'lo'`` can be chosen.
   Default is ``'tr'``.

- **ae_shape1/ae_shape2** (character, 0d/3d)
   Shape of the first/second pulse.

  - ``'impulse'``: Impulsive fields.
  - ``'Acos2'``: Envelope of cos\ :sup:`2`\ for a vector potential.
  - ``'Ecos2'``: Envelope of cos\ :sup:`2`\ for a scalar potential.

    If ``&system/iperiodic`` is ``3``, following can be also chosen,

  - ``'Acos3'``, ``'Acos4'``, ``'Acos6'``, and ``'Acos8'``: Envelopes of cos\ :sup:`3`\,cos\ :sup:`4`\, cos\ :sup:`6`\, and cos\ :sup:`8`\ for vector potentials.
  - (Trial) ``'Esin2sin'``, ``'Asin2cos'``, ``'Asin2cw'``, ``'input'``, and ``'none'`` can be also chosen.


- **e_impulse** (real(8), 0d/3d)
   Momentum of impulsive perturbation.
   This valiable has the dimention of momentum, energy*time/length.
   Default value is ``1d-2`` a.u.

..
 - **t_impulse**
   not yet implemented XXXX
..

   
- **E_amplitude1/E_amplitude2** (real(8), 0d/3d)
   Maximum amplitude of electric fields for the first/second pulse.
   This valiable has the dimension of electric field, energy/(length*charge).
   This valiable cannot be set with ``&emfield/rlaser_int_wcm2_1`` (``rlaser_int_wcm2_2``) simultaneously.

- **I_wcm2_1/I_wcm2_2** (real(8), 0d/3d)
   Peak laser intensity (W/cm\ :sup:`2`\) of the first/second pulse.
   This valiable cannot be set with ``&emfield/amplitude1`` (``amplitude2``) simultaneously.

- **tw1/tw2** (real(8), 0d/3d)
   Duration of the first/second pulse. Unit of time can be chosend 
   by ``&units/unit_time``.

- **omega1/omega2** (real(8), 0d/3d)
   Mean photon energy (average frequency multiplied by the Planck constant) of the first/second pulse. Unit of energy can be chosend 
   by ``&units/unit_energy``.

- **epdir_re1(3)/epdir_re2(3)** (real(8), 0d/3d)
   Real part of polarization vector for the first/second pulse.

- **epdir_im1(3)/epdir_im2(3)** (real(8), 0d/3d)
   Imaginary part of polarization vector for the first/second pulse.

- **phi_cep1/phi_cep2** (real(8), 0d/3d)
   Carrier emvelope phase of the first/second pulse.
   Default is ``0d0/0d0``.

- **t1_start** (real(8), 3d)
   Time-delay of the first pulse.
   Unit of time can be chosen by ``&units/unit_time``.
   (this is not available for multiscale option).
   Default is ``0d0``.

- **t1_t2** (real(8), 0d/3d)
   Time-delay between the first and the second pulses.
   Unit of time can be chosen by ``&units/unit_time``.

- (Trial) **yn_local_field** (character, 0d)
   The pulse is applied to a specific domain.
   Default is ``'n'``.

- **num_dipole_source** (integer, 0d)
   Number of radiation sources for optical near fields.
   Maximum number is ``2``.

- **vec_dipole_source(3,num_dipole_source)** (real(8), 0d)
   Dipole vectors of the radiation sources for the optical near fields.
   Unit of length can be chosen by ``&units/unit_length``.

- **cood_dipole_source(3,num_dipole_source)** (real(8), 0d)
   Central coordinates of the dipole vectors of the radiation sources.
   Unit of length can be chosen by ``&units/unit_length``.

- **rad_dipole_diele** (real(8), 0d)
   Radii of dielectric spheres for the radiation sources.
   Unit of length can be chosen by ``&units/unit_length``.



&multiscale
-----------

- (Trial) **fdtddim** (character, 3d)
   Dimension of FDTD calculation for multi-scale Maxwell-Kohn-Sham method.
   Default value is ``'1D'``. 

- (Trial) **twod_shape** (character, 3d)
   Boundary condision of the second dimension for FDTD calculation with 
   multi-scale Maxwell-Kohn-Sham method.
   Default value is ``'periodic'``.

- **nx_m** (integer, 3d)
   Number of macroscopic grid points inside materials for x-direction.

- (Trial) **ny_m/nz_m** (integer, 3d)
   Number of macroscopic grid points inside materials for (y/z)-direction.

- **hx_m** (real(8), 3d)
   Spacing of macroscopic grid points inside materials for (x)-direction.
   Unit of length can be chosen by ``&units/unit_length``.

- (Trial) **hy_m/hz_m** (real(8), 3d)
   Spacing of macroscopic grid points inside materials for (y/z)-direction.
   Unit of length can be chosen by ``&units/unit_length``.

- **nxvacl_m/nxvacr_m** (integer, 3d)
   Number of macroscopic grid points for vacumm region.
   ``nxvacl_m`` gives the number for negative x-direction in front of material,
   while ``nxvacr_m`` gives the number for positive x-direction behind the material.

- (Trial) **nx_origin_m/ny_origin_m/nz_origin_m** (integer, 3d)
   Origin coordinat of the grid points.
   Default value is ``'1'``.

- (Trial) **set_ini_coor_vel** (character, 3d)
   Set initial atomic coordinates and velocities for each macro-grid point. This must be given with specific directories and files: 
   Prepare ``directory``/multiscale/MXXXXXX/ini_coor_vel.dat, where 'XXXXXX' is the index number of the macro-grid point of the material region usually starting from '000001' up to the number of macro-grid point. The format of the file 'ini_coor_vel.dat' is just Rx, Ry, Rz, Vx, Vy, Vz (with space separation) for each atom (i.e. for each line), where the unit of the coordinates, Rx, Ry, Rz, is angstrom or a.u. speficied by ``unit_system`` but that of velocities is always a.u.. This option should be used together with ``read_gs_wfn_k_ms`` which is the option to read the ground state wave function for each macro-grid point. 
   Default value is ``'n'``.

- (Trial) **nmacro_write_group** (integer, 3d)
   If the number of macroscopic grids are very large, computers can be unstable by writing all information of all macroscopic grid points at the same time. To avoid that, the writings are divided by specifying this option. Writings will be done by each ``nmacro_write_group`` macroscopic grid points. (this number must be aliquot part of the total number of macroscopic grid points)
   Default value is ``'-1'``.

- (Trial) **file_macropoint** (character, 3d)
   If file name is specified in the option, the coordinates of the macropoints are set from the file.
   Default value is ``''``.


&maxwell
--------

- **al_em(3)** (real(8), 0d/3d)
   Size of simulation box in electromagnetic analysis. Unit of the length can be chosen by ``&units/unit_system``.

- **dl_em(3)** (real(8), 0d/3d)
   Spacing of real-space grids in electromagnetic analysis. Unit of length can be chosen by ``&units/unit_system``.

- **dt_em** (real(8), 0d/3d)
   Time step in electromagnetic analysis. Unit of time can be chosen by ``&units/unit_system``.

- **nt_em** (integer, 0d/3d)
   Number of total time steps for real-time propagation in electromagnetic analysis.

- **boundary_em(3,2)** (character, 0d/3d)
   Boundary condition in electromagnetic analysis. The first index(1-3 rows) corresponds to x, y, and z axes. The second index(1-2 columns) corresponds to bottom and top of the axes.  Default is ``'default'``. If ``&system/iperiodic=0``, ``'default'``, ``'pml'``, and ``'pec'`` can be chosen. ``'pml'`` is absorbing boundary and ``'pec'`` is perfect electric conductor. ``'default'`` is ``'pml'``. If ``&system/iperiodic=3``, ``'default'``, ``'pml'``, and ``'periodic'`` can be chosen. ``'periodic'`` is periodic boundary. ``'default'`` is ``'periodic'``.

- **shape_file** (character, 0d/3d)
   Name of shape file in electromagnetic analysis. The shape files can be generated by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **imedia_num** (integer, 0d/3d)
   Number of media in electromagnetic analysis. Default is ``0``.

- **type_media(:)** (character, 0d/3d)
   Type of media in electromagnetic analysis. ``'vacuum'``, ``'constant media'``, ``'pec'``, and ``'lorentz-drude'`` can be chosen. Default is ``'vacuum'``. If ``'lorentz-drude'`` is chosen, linear response calculation can be done by ``&emfield/ae_shape1 or ae_shape2='impulse'``.

- **epsilon(:)** (real(8), 0d/3d)
   Relative permittivity of the media in electromagnetic analysis. Default is ``1d0``.

- **rmu(:)** (real(8), 0d/3d)
   Relative permeability of the media in electromagnetic analysis. Default is ``1d0``.

- **sigma(:)** (real(8), 0d/3d)
   Conductivity of the media in electromagnetic analysis. Default is ``0d0``.

- **pole_num_ld(:)** (integer, 0d/3d)
   Number of poles of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. Default is ``1``.

- **omega_p_ld(:)** (real(8), 0d/3d)
   Plasma frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. Default is ``0d0``.

- **f_ld(:,:)** (real(8), 0d/3d)
   Oscillator strength of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``imedia_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``. Default is ``0d0``.

- **gamma_ld(:,:)** (real(8), 0d/3d)
   Collision frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``imedia_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``. Default is ``0d0``.

- **omega_ld(:,:)** (real(8), 0d/3d)
   Oscillator frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``imedia_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``. Default is ``0d0``.

- **wave_input** (character, 0d/3d)
   If ``'source'``, the incident pulse in electromagnetic analysis is generated by the incident current source. Default is ``'none'``.

- **ek_dir1(3)/ek_dir2(3)** (real(8), 0d/3d)
   Propagation direction of the first/second pulse.

- **source_loc1(3)/source_loc2(3)** (real(8), 0d/3d)
   Location of the incident current source of the first/second pulse. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/iperiodic=0`` while ranges from ``0`` to ``al_em`` for ``&system/iperiodic=3``.

- **iobs_num_em** (integer, 0d/3d)
   Number of observation point in electromagnetic analysis. Default is ``0``. From the obtained results, figure and animation files can be generated by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **iobs_samp_em** (integer, 0d/3d)
   Sampling time-step of the observation in electromagnetic analysis. Default is ``1``.

- **obs_loc_em(:,3)** (real(8), 0d/3d)
   Location of the observation point in electromagnetic analysis. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/iperiodic=0`` while ranges from ``0`` to ``al_em`` for ``&system/iperiodic=3``.

- **obs_plane_em(:)** (character, 0d/3d)
   Enable(``'y'``)/disable(``'n'``). Output of the electrmagnetic fields on the planes (xy, yz, and xz planes) for each observation point. This option must be ``'y'`` for generating animation files by using SALMON utilities (https://salmon-tddft.jp/utilities.html). Default is ``'n'``.

- (Trial) **wf_em** (character, 0d/3d)
   Enable(``'y'``)/disable(``'n'``). Applying a window function for linear response calculation when ``&calculation/theory=Maxwell``. Default is ``'y'``.

&analysis
---------

- **projection_option** (character, 3d)
   Methods of projection.
   
  - ``'no'``: no projection.
  - ``'gs'``: projection to eigenstates of ground-state Hamiltonian.
  - ``'rt'``: projection to eigenstates of instantaneous Hamiltonian.
  

- (Trial) **projection_decomp** (character, 3d)
   If ``'atom'`` combined with ``projection_option='gs'``, 
   the number of excited electron is decomposed into each atom 
   (this is printed in ``SYSname``\_nex_atom.data).
   Default is ``'n'``.

- **out_projection_step** (integer, 3d)
   Interval time step of projection analysis 
   if ``projection_option`` is not ``'no'``.
   Default is ``100``.

- **nenergy** (integer, 0d/3d)
   Number of energy grids for frequency-domain analysis.
   This parameter is required when `'impulse'` is choosen in `&emfield/ae_shape1|2`.

- **de** (real(8), 0d/3d)
   Energy spacing for analysis.
   Unit of energy can be chosen by ``&units/unit_energy``
   This parameter is required when `'impulse'` is choosen in `&emfield/ae_shape1|2`.

- **yn_out_psi** (character, 0d/3d)
   If ``'y'``, wavefunctions are output.
   For periodic system (``iperiodic=3``), it works only for ground state calculation. The converged wave functions of all orbitals with all k-points are printed in gs_wfn_cube or gs_wfn_vtk directory. The format is speficied by ``format3d``. 
   Default is ``'n'``.

- **yn_out_dos** (character, 0d/3d)
   If ``'y'``, density of state is output.
   Default is ``'n'``.

- **out_dos_start** (real(8), 0d/3d)
   Lower bound (energy) of the density of state spectra.
   If this value is lower than a specific value near the lowest energy level, 
   this value is overwritten by that value. 
   Default value is ``-1.d10`` eV.

- **out_dos_end** (real(8), 0d/3d)
   Upper bound (energy) of the density of state spectra.
   If this value is higher than a specific value near the highest energy level, 
   this value is overwritten by that value. 
   Default value is ``1.d10`` eV.

- **out_dos_nenergy** (integer, 0d/3d)
   Number of  energy points sampled in the density of state spectra.
   Default is ``601``.

- **out_dos_width** (real(8), 0d/3d)
   Smearing width used in the density of state spectra..
   Default is ``0.1`` eV.

- **out_dos_function** (character, 0d/3d)
   Choise of smearing method for the density of state spectra..
   ``gaussian`` and ``lorentzian`` function are available.
   Default is ``gaussian``.

- **yn_out_dos_set_fe_origin** (character, 0d/3d)
   If ``'y'``, the electron energy is shifted to fix the Fermi energy as zero point.
   For ``&system/iperiodic`` is ``0``, `` out_dos_fshift`` is not used 
   if ``&system/nstate`` is equal to ``&system/nelec``/2.
   Default is ``'n'``.

- **yn_out_pdos** (character, 0d)
   If ``'y'``, projected density of state is output.
   Default is ``'n'``.

- **yn_out_dns** (character, 0d/3d)
   If ``'y'``, the spatial electron density distribution at the ground state is output.
   Default is ``'n'``.

- **yn_out_dns_rt/out_dns_rt_step** ``Character/Integer``; 0d/3d)
   If ``'y'``,  the spatiotemporal electron density distribution during real-time time-propagation is output
   every ``outdns_rt_step`` time steps.
   Default is ``'n'``.

- (Trial) **yn_out_dns_trans/out_dns_trans_energy** ``Character/Real(8)``; 3d)
   If ``'y'``, transition in different density from the ground state at specified field frequency omega(given by ``out_dns_trans_energy``) is calculated by drho(r,omega)=FT(rho(r,t)-rho_gs(r))/T.
   Default is ``'n'/1.55eV``.

- **yn_out_elf** (character, 0d)
   If ``'y'``, electron localization function is output.
   Default is ``'n'``.

- **yn_out_elf_rt/out_elf_rt_step** ``Character/Integer``; 0d)
   If ``'y'``, electron localization function 
   during real-time time-propagation is output
   every ``out_elf_rt_step`` time steps.
   Default is ``'n'``.

- **yn_out_estatic_rt/out_estatic_rt_step** ``Character/Integer``; 0d)
   If ``'y'``, static electric field
   during real-time time-propagation is output
   every ``out_estatic_rt_step`` time steps.
   Default is ``'n'``.

- (Trial) **yn_out_rvf_rt/out_rvf_rt_step** ``Character/Integer``; 3d)
   If ``'y'``, coordinates[A], velocities[au], forces[au] on atoms
   during real-time time-propagation are printed in ``SYSname``\_trj.xyz
   every ``out_rvf_rt_step`` time steps.
   If ``use_ehrenfest_md='y'``, 
   the printing option is automatically turned on.
   Defaults are ``'n'/10``.

- (Trial) **yn_out_tm** (character, 3d)
   If ``'y'``, transition moments between occupied and virtual orbitals are printed into ``SYSname``\_tm.data after the ground state calculation.
   Defaults are ``'n'``.

- **format_voxel_data** (character, 0d/3d)
   File format for three-dimensional volumetric data.
   ``'avs'``, ``'cube'``, and ``'vtk'`` can be chosen.
   Default is ``'cube'``.

- **nsplit_voxel_data** (integer, 0d)
   Number of separated files for three dimensional data.
   Effective only when ``format3d`` is ``'avs'``.
   ``numfiles_out_3d`` must be less than or equal to number of processes.
   Default is ``1``.

- (Trial) **timer_process** (character, 0d)
   Basically, elapsed times are written in the output file. 
   But if ``timer_process`` is ``'y'``, 
   files of elapsed times for every process are also generated. 
   This variable is effective only for the real-time caululation.
   Default is ``'n'``.


&poisson
--------

- **layout_multipole** (character, 0d)
   A variable to determine how to put multipoles in the Hartree potential calculation. Default is ``3``.

  - ``1``: A single pole is put at the center.
  - ``2``: Multipoles are put at the center of atoms.
  - ``3``: Multipoles are put at the center of mass of electrons in prepared cuboids.

- **num_multipole_xyz(3)** (integer, 0d)
   Number of multipoles when ``meo`` is ``3``. Default is ``0,0,0``. When default is set, number of multipoles is calculated automatically.


&ewald
------

- **newald** (integer, 3d)
   Parameter for Ewald method. 
   Short-range part of Ewald sum is calculated within ``newald`` th
   nearlist neighbor cells.
   Default is ``4``.

- **aewald** (real(8), 3d)
   Square of range separation parameter for Ewald method in atomic unit. 
   Default is ``0.5``.



&opt(Trial)
-------------

- **nopt** (integer, 0d/3d)
   xxx

- (Trial) **convrg_opt_fmax** (real(8), 3d)
   Convergence threshold of optimization in maximum force.
   Default is ``1d-3``.

..  
  - (Trial) **cg_alpha_up** (real(8), 3d)
    Parameter for up-rate of step length in line search in conjugated gradient method.
    Default is ``1.3``.

  - (Trial) **cg_alpha_down** (real(8), 3d)
    Parameter for down-rate of step length in line search in conjugated gradient method.
    Default is ``0.5``.

  - (Trial) **cg_alpha_ini** (real(8), 3d)
    Parameter for initial step length in line search in conjugated gradient method. (currently not available)
    Default is ``0.8``.

  - (Trial) **convrg_scf_ene** (real(8), 3d)
    Convergence threshold of ground state SCF calculation in energy difference at each optimization step. If negative number no threshold (SCF loop is up to ``Nscf``). The other SCF thresholds such as ``threshold`` in ``&scf`` are also applied (if you do not want to use it, set very small number). 
    Default is ``-1.0``.

  - (Trial) **convrg_scf_force** (real(8), 3d)
    Convergence threshold of ground state SCF calculation in force (average over atoms) difference. If negative number no threshold (SCF loop is up to ``Nscf``). The other SCF thresholds such as ``threshold`` in ``&scf`` are also applied (if you do not want to use it, set very small number). 
    Default is ``-1.0``.

  - (Trial) **convrg_opt_ene** (real(8), 3d)
    Convergence threshold of optimization in energy difference. (currently not available)
    Default is ``1d-6``.
..


&md(Trial)
-----------
- (Trial) **ensemble** (character, 3d)
   Ensemble in MD option: "NVE" or "NVT".
   Default is ``"NVE"``.

- (Trial) **thermostat** (character, 3d)
   Thermostat in "NVT" option: (currently only ``nose-hoover``).
   Default is ``"nose-hoover"``.

- (Trial) **step_velocity_scaling** (integer, 3d)
   Time step interval for velocity-scaling. Velocity-scaling is applied if this is set to positive.
   Default is ``-1``.

- (Trial) **step_update_ps/step_update_ps2** ``Integer/Integer``; 3d)
   Time step interval for updating pseudopotential (Larger number makes calculation time reduce greatly, but gets inaccurate) in case of ``use_ehrenfest_md=y``. ``step_update_ps`` is for full update and ``step_update_ps2`` is for update without changing grid points array.
   Default is ``10/1``.

- (Trial) **temperature0_ion_k** (real(8), 3d)
   Setting temperature [K] for NVT ensemble, velocity scaling and generating initial velocities.
   Default is ``298.15``.

- (Trial) **yn_set_ini_velocity** (character, 3d)
   Initial velocities are set.
   Default is ``n``.

  - ``y``: Generate initial velocity with Maxwell-Bortzman distribution.
  - ``r``: Read initial velocity from file specified by keyword of ``file_ini_velocity``. This is, for example, used for restarting MD from the previous run. The last atomic coordinates and velocities are printed in ``SYSname``\_trj.xyz. (atomic coordinate also should be copied from the previous output and put in the next input file for restart)

    
- (Trial) **file_ini_velocity** (character, 3d)
   File name for initial velocities. This is read when ``set_ini_velocity`` is ``'r'``. The format is simply vx(iatom) vy(iatom) vz(iatom) in each line. The order of atoms must be the same as the given coordinates in the main input file. In case of using nose-hoover thermostat, a thermostat variable should be put at the last line (all atomic unit). 
   Default is ``none``.

- (Trial) **seed_ini_velocity** (integer, 3d)
   Random seed (integer number) to generate initial velocity if ``set_ini_velocity`` is set to y.
   Default is ``123``.

- (Trial) **thermostat_tau** (real(8), 3d)
   Parameter in Nose-Hoover method: controlling time constant for temperature.
   Default is ``41.34[au] or 1.0[fs]``.

- (Trial) **yn_stop_system_momt** (character, 3d)
   Center of mass is stopped every time step.
   Default is ``n``.


&code
-----

- **yn_want_stencil_openmp_parallelization(yn)**

- **yn_want_stencil_hand_vectorization(yn)**

- **yn_force_stencil_openmp_parallelization(yn)**

- **yn_force_stencil_sequential_computation(yn)**

- **yn_want_communication_overlapping(yn)**

   

**Following variables are moved from the isolated part. Some of them may be added to common input, be combined to it, and be removed.**


&group_fundamental(Trial)
-------------------------

- (Trial) **iditer_nosubspace_diag** (integer, 0d)
   Iterations for which subspace diagonalization is not done if ``&scf/subspace_diagonalization`` is ``'y'``.
   Default is ``10``.

- (Trial) **ntmg** (integer, 0d)
   Number of multigrid calculation for gs. At the moment, there is a malfunction in this variable, and recovery is needed.
   Default is ``1``.

- (Trial) **idisnum(2)** (integer, 0d)
   Label numbers for two atoms which are measured the distance. 
   Default is ``(/1,2/)``.

- (Trial) **iwrite_projection** (integer, 0d)
   A variable for projection. 
   Default is ``0``.

- (Trial) **itwproj** (integer, 0d)
   The projection is calculated every ``itwproj`` time steps. 
   Default is ``-1``.

- (Trial) **iwrite_projnum** (integer, 0d)
   There is a malfunction in this variable.

- (Trial) **itcalc_ene** (integer, 0d)
   Total energy is calculated every ``itcalc_ene`` time steps. There may be a malfunction in this variable.
   Default is ``1``.


&group_parallel(Trial)
-----------------------

- (Trial) **isequential** (integer, 0d)
   A variable to determine the way of assignment of processes.
   Default is ``2``.

- (Trial) **imesh_s_all** (integer, 0d)
   A variable to determine how to use processes if total number of processes 
   and number of processes for Hartree/Exc calculation differ. 
   There may be a malfunction in this variable.
   Default is ``1``.

- (Trial) **iflag_comm_rho** (integer, 0d)
   This variable may be removed. 


&group_hartree(Trial)
----------------------

- (Trial) **hconv** (real(8), 0d)
   A convergence value for the Hartree-cg calculation. 
   The convergence is checked by ||tVh(i)-tVh(i-1)||\ :sup:`2`\/(number of grids).
   Default is ``1d-15`` a.u. (= 1.10d-13 A\ :sup:`3`\eV\ :sup:`2`\)

- (Trial) **lmax_meo** (integer, 0d)
   A maximum angular momentum for multipole expansion in the Hartree-cg calculation. 
   Default is ``4``.



&group_file(Trial)
-------------------

- (Trial) **ic** (integer, 0d)
   A variable to check whether reentrance is done or not in the ground state calculation. 
   Default is ``0``.

- (Trial) **oc** (integer, 0d)
   A variable to check whether intermediate files are generated in the ground state calculation. 
   Default is ``1``.

- (Trial) **ic_rt** (integer, 0d)
   A variable to check whether reentrance is done or not in the time propagation calculation. 
   Default is ``0``.

- (Trial) **oc_rt** (integer, 0d)
   A variable to check whether intermediate files are generated in the time propagation calculation. 
   Default is ``0``.


&group_others(Trial)
---------------------

- (Trial) **iparaway_ob** (integer, 0d)
   A variable to determine the way of division for orbitals. 
   ``1`` is block division, and ``2`` is cyclic division.
   Default is ``2``.

- (Trial) **iswitch_orbital_mesh** (integer, 0d)
   A variable to apply descending order for orbitals in the ground state calculation.
   Default is ``0``.

- (Trial) **iflag_psicube** (integer, 0d)
   A variable to generate cube files for wave functions. This variable will be removed.

- (Trial) **file_ini** (character, 0d)
   A input file to align wavefunctions. 
   Default is ``'file_ini'``.

- (Trial) **num_projection** ``Interger``; 0d)
   Number of orbitals for projections.
   Default is ``1``.

- (Trial) **iwrite_projection_ob(200)** ``Interger``; 0d)
   Orbital number to be written as projections.
   Default is ``(1/2/3/.../200)``.

- (Trial) **iwrite_projection_k(200)** ``Interger``; 0d)
   This variable will be removed.

- (Trial) **filename_pot** (character, 0d)
   Name of file to be written local potentials. 
   Default is ``'pot'``.

- (Trial) **iwrite_external** (integer, 0d)
   A variable to generate file to be written local potentials. 
   Default is ``0``.

- (Trial) **iflag_dip2** (integer, 0d)
   A variable to determine whether dipole moments are calculated in divided area. 
   Default is ``0``.

- (Trial) **iflag_intelectron** (integer, 0d)
   A variable related to the quadrupole caluclation.
   Default is ``0``.

- (Trial) **num_dip2** (integer, 0d)
   Number of area where dipole moments are calculated.
   Default is ``1``.

- (Trial) **dip2boundary(100)** (real(8), 0d)
   Boundary position of area where dipole moments are calculated.
   Default is ``0`` a.u.

- (Trial) **dip2center(100)** (real(8), 0d)
   Origin in the dipole moment calculation. 
   Default is ``0`` a.u.

- (Trial) **iflag_fourier_omega** ``integer``; 0d)
   A variable to determine whether Fourier transformation of 3d data for difference of density is calclated. 
   Default is ``0``.

- (Trial) **num_fourier_omega** (integer, 0d)
   Number of energies for which the Fourier transformation is calclated. 
   Default is ``1``.

- (Trial) **fourier_omega(200)** (real(8), 0d)
   Energies for which the Fourier transformation is calclated. 
   Default is ``0`` a.u.

- (Trial) **itotntime2** (integer, 0d)
   Number of time steps in the reentrance for real-time calculation.
   There may be a malfunction in this variable.
   Default is ``0``.

- (Trial) **iwdenoption** (integer, 0d)
   A variable to determine whether 3d output is generated in real-time calculation. 
   This variable will be removed.

- (Trial) **iwdenstep** (integer, 0d)
   3d output is generated every ``iwdenstep`` time steps.
   This variable will be removed.

- (Trial) **iflag_estatic** (integer, 0d)
   A variable to determine whether 3d output for the static electric field is generated in real-time calculation. 
   This variable will be removed.


   
.. _&calculation: #calculation
.. _&control: #control
.. _&units: #units
.. _&parallel: #parallel
.. _&system: #system
.. _&atomic-red-coor: #atomic_red_coor
.. _&atomic-coor: #atomic_coor
.. _&pseudo: #pseudo
.. _&functional: #functional
.. _&rgrid: #rgrid
.. _&kgrid: #kgrid
.. _&tgrid: #tgrid
.. _&propagation: #propagation
.. _&scf: #scf
.. _&emfield: #emfield
.. _&multiscale: #multiscale
.. _&maxwell: #maxwell
.. _&analysis: #analysis
.. _&poisson: #poisson
.. _&ewald: #ewald
.. _&opt: #opt
.. _&md: #md
.. _&group_fundamental: #group_fundamental
.. _&group_parallel: #group_parallel
.. _&group_hartree: #group_hartree
.. _&group_file: #group_file
.. _&group_others: #group_others



