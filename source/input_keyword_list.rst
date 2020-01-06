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
-  `&opt[Trial]`_
-  `&md[Trial]`_
-  `&group_fundamental[Trial]`_
-  `&group_parallel[Trial]`_  
-  `&group_hartree[Trial]`_ 
-  `&group_file[Trial]`_
-  `&group_others[Trial]`_


&calculation
------------

- **theory** (character, default='none')
   | Choice of Calculation theories.
   |  Options
   |    ``DFT``  / ground state calculation based on DFT
   |    ``DFT_MD``  / adiabatic ab initio MD simulations based on DFT
   |    ``TDDFT_response``  / simulations based on TDDFT for response
   |    ``TDDFT_pulse``  / simulations based on TDDFT using pulsed light
   |    ``Single_scale_Maxwell_TDDFT``  / coupled Maxwell and TDDFT single-scale simulation
   |    ``Multi_scale_Maxwell_TDDFT``  / coupled Maxwell and TDDFT multi-scale simulation
   |    ``Maxwell``  / electromagnetic simulations based on the Maxwell's equations (o

- **yn_md** (character, Default='n')[Trial]
   | Available for ``theory='DFT'`` (Adiabatic ground-state MD) and ``theory='TDDFT_pulse'`` (Ehrenfest MD).
   Enable(``'y'``)/disable(``'n'``). Molecular dynamics option.

- **yn_opt** (character, Default='n')[Trial]
   | Available for ``theory='DFT'``.
   Enable(``'y'``)/disable(``'n'``). Geometry optimization option.


&control
--------

- **sysname** (character, Default='default')
   | Available for ``theory='xxx'``.
   Name of calculation. This is used for a prefix of output files.

- **base_directory** (character, Default='./')
   | Available for ``theory='xxx'``.
   Name of a default directory, where the basic results will be written down.

- **output_buffer_interval** (integer, Default=-1)
   | Available for ``theory='xxx'``.
   xxx.

- **yn_restart** (character, Default='new')[Trial]
   | Available for ``theory='xxx'``.
   Flag for restart, ``'new'`` or ``'restart'``.

- **directory_read_data** (character, Default='restart/')[Trial]
   | Available for ``theory='xxx'``.
   xxx.

- **yn_self_checkpoint** (character, Default='n')[Trial]
   | Available for ``theory='xxx'``.
   xxx.

- **checkpoint_interval** (integer, Default=0)[Trial]
   | Available for ``theory='xxx'``.
   Frequency of backup during the time-propagation. 
   If ``0`` is set, the backup is not performed.

- **time_shutdown** (real(8), Default=-1d0)[Trial]
   | Available for ``theory='xxx'``.
   Timer for automatic shutdown. The unit is always second.
   If negative time is chosen, the automatic shutdown is not performed.

&units
------

- **unit_system** (character, Default='au')
   Unit for input variables. 
   If ``'au'`` or ``'a.u.'``, atomic unit system is used. 
   If ``'A_eV_fs'``, Angstrom-eV-fs unit system is used. 


&parallel
---------

- **yn_domain_parallel** (character, Default='n')[Trial]
   | Old infomation: 3d
   If specified ``yn_domain_parallel='y'`` and ``&system/yn_periodic='y'``, program codes for domain parallel version run in periodic system calculations.

- **nproc_k/nproc_ob/nproc_domain_orbital(3)/nproc_domain_general(3)** (integer, Default=0)
   | Old infomation: 0d
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

- **yn_ffte** (character, Default='ft')
   | Old infomation: 0d
   Method of Fourier transformation.  ``'ft'``,  ``'FT'``, ``'ffte'`` or ``'FFTE'`` can be chosen.
   This variable is effective only when ``yn_domain_parallel='y'`` and ``&system/yn_periodic='y'``.

- **process_allocation** (character, Default='grid_sequential')
   | Old infomation: 0d
   xxx.


&system 
-------

- **yn_periodic** (character, Default='n')
   Dimension for periodic boundary condition.
   ``'n'`` is for isolated systems, and 
   ``'y'`` is for solids.

- **ispin** (integer, Default=0)
   | Old infomation: 0d
   Variable for classification of closed shell systems and open shell systems.
   ``0`` is for closed shell systems, and
   ``1`` is for open shell systems.

- **al(3)** (real(8), Default=0d0)
   Lattice constants. Unit of the length can be chosen by ``&units/unit_system``.

- **al_vec1(3)/al_vec2(3)/al_vec3(3)** (real(8), Default=0d0)
   | Old infomation: 3d
   xxx.

- **isym** (integer, Default=0)
   | Old infomation: 3d
   Number of symmetries that can be used for reduction of k-points.

- **crystal_structure** (character, Default='none')
   | Old infomation: 3d
   Name of symmetry that can be used for the reduction of # of k-points.

- **nstate** (integer, Default=0)
   Number of states/bands.

- **nstate_spin(2)** (integer, Default=0)
   | Old infomation: 0d
   Number of states/bands can be specified independently by ``nstate_spin(1)/nstate_spin(2)``.
   This option is incompatible with ``nstate``

- **nelec** (integer, Default=0)
   Number of valence electrons.

- **nelec_spin(2)** (integer, Default=0)
   | Old infomation: 0d
   Number of up/down-spin electrons can be specified independently by ``nelec_spin(1)/nelec_spin(2)``.
   This option is incompatible with ``nelec``

- **temperature** (real(8), Default=-1d0)
   | Old infomation: 3d
   Temperature of electrons. When you calculate a system of zero band-gap energy like metals, zero or positive number of the temperature should be given.
   Unit of the energy can be chosen ``&units/unit_system``. 
   Default is ``-1.0`` (this is for system which has a band gap energy).

- **temperature_k** (real(8), Default=-1d0)[Trial]
   | Old infomation: 0d
   Temperature of electrons [K]. Default is ``-1.d0``.

- **nelem** (integer, Default=0)
   Number of elements that will be used in calculations.

- **natom** (integer, Default=0)
   Number of atoms in a calculation cell.


- **file_atom_red_coor** (character, Default='none')[Trial]
   | Old infomation: 3d
   File name of atomic positions. In this file, 
   the atomic coordinates can be written in reduced coordinates.
   This option is incompatible with 
   ``&system/file_atom_coor``,
   ``&atomic_coor``, and 
   ``&atomic_red_coor``.

- **file_atom_coor** (character, Default='none')[Trial]
   | Old infomation: 0d
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

- **file_pseudo(:)** (character, Default='none')
   Name of pseudopotential files.

- **lmax_ps(:)** (integer, Default=-1)
   Maximum angular momentum of pseudopotential projectors.

- **lloc_ps(:)** (integer, Default=-1)
   Angular momentum of pseudopotential that will be treated as local.

- **izatom(:)** (integer, Default=-1)
   Atomic number.

- **yn_psmask(:)** (character, Default='n')[Trial]
   Enable(``'y'``)/disable(``'n'``) 
   Fourier filtering for pseudopotentials. 

- **alpha_mask(:)** (real(8), Default=0.8d0)[Trial]
   Parameter for the Fourier filtering for pseudopotential.

- **gamma_mask(:)** (real(8), Default=1.8d0)[Trial]
   Parameter for the Fourier filtering for pseudopotential.

- **eta_mask(:)** (real(8), Default=15.0d0)[Trial]
   Parameter for the Fourier filtering for pseudopotential.


&functional
-----------

- **xc** (character, Default='none')
   Exchange-correlation functionals.
   At present version, the functional 'PZ', 'PZM' and 'TBmBJ' is available for both 0d/3d calculations, and the functionals 'TPSS' and 'VS98' are available for 3d calculations.

  - ``'PZ'``: Perdew-Zunger LDA :Phys. Rev. B 23, 5048 (1981).
  - ``'PZM'``: Perdew-Zunger LDA with modification to improve sooth connection between high density form and low density one. :J. P. Perdew and Alex Zunger, Phys. Rev. B 23, 5048 (1981).
  - ``'TBmBJ'``: Tran-Blaha meta-GGA exchange with Perdew-Wang correlation. :Fabien Tran and Peter Blaha, Phys. Rev. Lett. 102, 226401 (2008). John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992).
  - ``'TPSS'``: Tao, Perdew, Staroverov and Scuseria meta-GGA exchange correlation. :J. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria, Phys. Rev. Lett. 91, 146401 (2003).
  - ``'VS98'``:  van Voorhis and Scuseria exchange with Perdew-Wang correlation: T. Van Voorhis and G. E. Scuseria, J. Chem. Phys. 109, 400 (1998).

- **cname, xname** (character, Default='none')
   xxx.

- **alibxc, alibx, alibc** (character, Default='none')
   By specifying ``alibxc``, the functionals prepared in libxc package are available. 
   They can be set indivisually by specifying ``alibx`` and ``alibc``.
   To use libxc libraries, ``--with-libxc`` option must be added in excecuting configure. 
   The available option of the exchange-correlation functionals are listed in the LibXC website. 
   [See http://www.tddft.org/programs/libxc/functionals/]
   
- **cval** (real(8), Default=-1d0)
   | Old infomation: 3d
   Mixing parameter in Tran-Blaha meta-GGA exchange potential. If ``cval`` is set to a minus value, the mixing-parameter computed
   by the formula in the original paper [Phys. Rev. Lett. 102, 226401 (2008)].
   Default is estimated from :math:`\left\langle |\nabla \rho(\mathbf{r};t)| / \rho(\mathbf{r};t) \right\rangle`.


&rgrid
------

- **dl(3)** (real(8), Default=0d0)
   Spacing of real-space grids. Unit of length can be chosen by
   ``&units/unit_system``.
   This valiable cannot be set with 
   ``&rgrid/num_rgrid`` simultaneously.
   If ``&system/yn_periodic`` is set to ``'y'``,
   the actual grid spacing is automatically refined in calculations
   so that the size of the simulation box
   ``&system/al(3)`` becomes divisible by the spacing.

- **num_rgrid(3)** (integer, Default=0)
   | Old infomation: 3d
   Number of real-space grids.
   This valiable cannot be set with 
   ``&rgrid/dl`` simultaneously.


&kgrid
------

- **num_kgrid(3)** (integer, Default=1)
   | Old infomation: 3d
   Number of k-points (grid points of k-vector) discretizing
   the Brillouin zone.
   Each component must be even.

- **file_kw** (character, Default='none')
   | Old infomation: 3d
   Name of a file for flexible k-point sampling.
   This file will be read if ``num_kgrid`` is smaller than 1.


&tgrid
------

- **nt** (integer, Default=0)
   Number of total time steps for real-time propagation.

- **dt** (real(8), Default=0d0)
   Time step. Unit of time can be chosen by ``&units/unit_system``.


&propagation
------------

- **n_hamil** (integer, Default=4)[Trial]
   | Old infomation: 0d
   Order of Taylor expansion of a propagation operator.

- **propagator** (character, Default=middlepoint')
   | Old infomation: 3d
   Choice of Propagator.
   ``middlepoint`` is an propagator
   with the Hamiltoinan at midpoint of two-times.
   ``etrs`` is enforced time-reversal symmetry propagator.
   [M.A.L. Marques, A. Castro, G.F. Bertsch, and A. Rubio, Comput. Phys. Commun., 151 60 (2003)].

- **yn_fix_func** (character(1), Default='n')[Trial]
   | Old infomation: 3d
   Option not to update functional (or Hamiltonian) in RT calculation, i.e., keep ground state Hamiltonian during time-evolution.

&scf
----

- **method_min** (character, Default='cg') 
   xxx.

- **ncg** (integer, Default=5)
   Number of interation of Conjugate-Gradient method for each scf-cycle.

- **method_mixing** (character, Default='broyden') 
   | Old infomation: 0d
   Methods for density/potential mixing for scf cycle. ``simple`` and ``broyden`` can be chosen.

- **mixrate** (real(8), Default=0.5d0)
   | Old infomation: 0d
  Mixing ratio for simple mixing.

- **nmemory_mb** (integer, Default=8)
   Number of stored densities at previous scf-cycles for 
   the modified-Broyden method. 
   If ``&system/yn_periodic`` is ``'n'``, ``nmemory_mb`` must be less than 21.

- **alpha_mb** (real(8), Default=0.75d0)
   Parameter of the modified-Broyden method.
 
- **fsset_option** (character, Default='n') 
   xxx.

- **nfsset_start** (integer, Default=75) 
   xxx.

- **nfsset_every** (integer, Default=25) 
   xxx.

- **nscf** (integer, Default=0)
   Number of maximum scf cycle.

- **yn_subspace_diagonalization** (character, Default='y')
   | Old infomation: 0d
   Enable(``'y'``)/disable(``'n'``) 
   subspace diagonalization during scf cycle.

- **convergence** (character, Default='rho_dne')
   Choice of quantity that is used for convergence check in a scf calculation. 

  - ``'rho_dne'``: Convergence is checked by sum_ix|rho(ix,iter)-rho(ix,iter-1)|dx/N, where iter is an iteration number of the scf calculation and N is ``&system/nelec``, the number of the valence electrons.

   For isolated systems, the followings can also be chosen.

  - ``'norm_rho'``: Convergence is checked by the square of the norm of difference of density, ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\=sum_ix|rho(ix,iter)-rho(ix,iter-1)|\ :sup:`2`\. 
  - ``'norm_rho_dng'``: Convergence is checked by ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\/(number of grids). "dng" means "devided by number of grids".
  - ``'norm_pot'``: Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\, where Vlocal is Vh + Vxc + Vps_local.
  - ``'pot_dng'``: Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\/(number of grids).

- **threshold** (real(8), Default=1d-17)
   Threshold for convergence check that is used when ``'rho_dne'`` is specified.
   Default is ``1d-17``. 
   XXX(threshold_norm_rho (real(8), Default=))XXX
   Threshold for convergence check that is used when either ``'norm_rho'`` or ``'norm_rho_dng'`` is specified. ``threshold_norm_rho`` must be set when either ``'norm_rho'`` or ``'norm_rho_dng'`` is specified.
   Default is ``-1d0`` a.u. (1 a.u.= 45.54 A\ :sup:`-6`\)
   XXX(threshold_norm_pot (real(8), Default=))XXX
   Threshold for convergence check that is used when either ``'norm_pot'`` or ``'norm_pot_dng'`` is specified. ``threshold_norm_pot`` must be set when either ``'norm_pot'`` or ``'norm_pot_dng'`` is specified.
   Default is ``-1d0`` a.u. (1 a.u.= 33.72x10\ :sup:`4`\ A\ :sup:`-6`\eV\ :sup:`2`\)

- **omp_loop** (character, Default='k')
   | Old infomation: 3d
   XXX only ARTED XXX
   Loop for OpenMP parallelization in the ground state SCF if periodic boundary system is used. 

  - ``k``: parallelization for k-point loop (Default).
  - ``b``: parallelization mainly for band orbital loop (sometimes space grid loop too). This works efficiently if the number of k-point treated in each node is small (e.x. the case of single k-point for each node)

- **skip_gsortho** (character, Default='n')[Trial]
   | Old infomation: 3d
   XXX only ARTED XXX
   Flag to skip Gram-Schmidt orthogonalization in CG loop if periodic boundary system is used. If this is skipped the more iteration number is necessary to get convergence but each iteration step gets faster. If ``omp_loop=b``, this flag is always applied.

- **iditer_notemperature** (integer, Default=10) 
   xxx.


&emfield
--------

- **trans_longi** (character, Default='tr')
   | Old infomation: 3d
   Geometry of solid-state calculations.
   Transverse ``'tr'`` and longitudinal ``'lo'`` can be chosen.

- **ae_shape1/ae_shape2** (character, Default='none')
   Shape of the first/second pulse.

  - ``'impulse'``: Impulsive fields.
  - ``'Acos2'``: Envelope of cos\ :sup:`2`\ for a vector potential.
  - ``'Ecos2'``: Envelope of cos\ :sup:`2`\ for a scalar potential.

    If ``&system/yn_periodic`` is ``'y'``, following can be also chosen,

  - ``'Acos3'``, ``'Acos4'``, ``'Acos6'``, and ``'Acos8'``: Envelopes of cos\ :sup:`3`\,cos\ :sup:`4`\, cos\ :sup:`6`\, and cos\ :sup:`8`\ for vector potentials.
  - [Trial] ``'Esin2sin'``, ``'Asin2cos'``, ``'Asin2cw'``, ``'input'``, and ``'none'`` can be also chosen.


- **e_impulse** (real(8), Default=1d-2 a.u.)
   Momentum of impulsive perturbation.
   This valiable has the dimention of momentum, energy*time/length.

..
 - **t_impulse**
   not yet implemented XXXX
..
   
- **E_amplitude1/E_amplitude2** (real(8), Default=0d0)
   Maximum amplitude of electric fields for the first/second pulse.
   This valiable has the dimension of electric field, energy/(length*charge).
   This valiable cannot be set with ``&emfield/I_wcm2_1`` (``I_wcm2_2``) simultaneously.

- **I_wcm2_1/I_wcm2_2** (real(8), Default=-1d0)
   Peak laser intensity (W/cm\ :sup:`2`\) of the first/second pulse.
   This valiable cannot be set with ``&emfield/E_amplitude1`` (``E_amplitude2``) simultaneously.

- **tw1/tw2** (real(8), Default=0d0)
   Duration of the first/second pulse. Unit of time can be chosend 
   by ``&units/unit_time``.

- **omega1/omega2** (real(8), Default=0d0)
   Mean photon energy (average frequency multiplied by the Planck constant) of the first/second pulse. Unit of energy can be chosend 
   by ``&units/unit_energy``.

- **epdir_re1(3)/epdir_re2(3)** (real(8), Default=1d0, 0d0, 0d0)
   Real part of polarization vector for the first/second pulse.

- **epdir_im1(3)/epdir_im2(3)** (real(8), Default=0d0)
   Imaginary part of polarization vector for the first/second pulse.

- **phi_cep1/phi_cep2** (real(8), Default=0d0)
   Carrier emvelope phase of the first/second pulse.
   Default is ``0d0/0d0``.

- **t1_t2** (real(8), Default=0d0)
   Time-delay between the first and the second pulses.
   Unit of time can be chosen by ``&units/unit_time``.

- **t1_start** (real(8), Default=0d0)
   | Old infomation: 3d
   Time-delay of the first pulse.
   Unit of time can be chosen by ``&units/unit_time``.
   (this is not available for multiscale option).

- **yn_local_field** (character, Default='n')[Trial]
   | Old infomation: 0d
   The pulse is applied to a specific domain.

- **rlaserbound_sta/rlaserbound_end** (real(8), Default=-1.d7 a.u.)
   xxx.

- **num_dipole_source** (integer, Default=0)
   | Old infomation: 0d
   Number of radiation sources for optical near fields.
   Maximum number is ``2``.

- **vec_dipole_source(3,num_dipole_source)** (real(8), Default=0d0)
   | Old infomation: 0d
   Dipole vectors of the radiation sources for the optical near fields.
   Unit of length can be chosen by ``&units/unit_length``.

- **cood_dipole_source(3,num_dipole_source)** (real(8), Default=0d0)
   | Old infomation: 0d
   Central coordinates of the dipole vectors of the radiation sources.
   Unit of length can be chosen by ``&units/unit_length``.

- **rad_dipole_diele** (real(8), Default=2d0 a.u.)
   | Old infomation: 0d
   Radii of dielectric spheres for the radiation sources.
   Unit of length can be chosen by ``&units/unit_length``.



&multiscale
-----------

- **fdtddim** (character, Default='1d')[Trial]
   | Old infomation: 3d
   Dimension of FDTD calculation for multi-scale Maxwell-Kohn-Sham method.

- **twod_shape** (character, Default='periodic')[Trial]
   | Old infomation: 3d
   Boundary condision of the second dimension for FDTD calculation with 
   multi-scale Maxwell-Kohn-Sham method.

- **nx_m** (integer, Default=1)
   | Old infomation: 3d
   Number of macroscopic grid points inside materials for x-direction.

- **ny_m/nz_m** (integer, Default=1)[Trial]
   | Old infomation: 3d
   Number of macroscopic grid points inside materials for (y/z)-direction.

- **hx_m** (real(8), Default=0d0)
   | Old infomation: 3d
   Spacing of macroscopic grid points inside materials for (x)-direction.
   Unit of length can be chosen by ``&units/unit_length``.

- **hy_m/hz_m** (real(8), Default=0d0)[Trial]
   | Old infomation: 3d
   Spacing of macroscopic grid points inside materials for (y/z)-direction.
   Unit of length can be chosen by ``&units/unit_length``.

- **nxvacl_m/nxvacr_m** (integer, Default=1/0)
   | Old infomation: 3d
   Number of macroscopic grid points for vacumm region.
   ``nxvacl_m`` gives the number for negative x-direction in front of material,
   while ``nxvacr_m`` gives the number for positive x-direction behind the material.

- **nx_origin_m/ny_origin_m/nz_origin_m** (integer, Default=1)[Trial]
   | Old infomation: 3d
   Origin coordinat of the grid points.

- **file_macropoint** (character, Default='')[Trial]
   | Old infomation: 3d
   If file name is specified in the option, the coordinates of the macropoints are set from the file.

- **set_ini_coor_vel** (character, Default='n')[Trial]
   | Old infomation: 3d
   Set initial atomic coordinates and velocities for each macro-grid point. This must be given with specific directories and files: 
   Prepare ``directory``/multiscale/MXXXXXX/ini_coor_vel.dat, where 'XXXXXX' is the index number of the macro-grid point of the material region usually starting from '000001' up to the number of macro-grid point. The format of the file 'ini_coor_vel.dat' is just Rx, Ry, Rz, Vx, Vy, Vz (with space separation) for each atom (i.e. for each line), where the unit of the coordinates, Rx, Ry, Rz, is angstrom or a.u. speficied by ``unit_system`` but that of velocities is always a.u.. This option should be used together with ``read_gs_wfn_k_ms`` which is the option to read the ground state wave function for each macro-grid point. 

- **nmacro_write_group** (integer, Default=-1)[Trial]
   | Old infomation: 3d
   If the number of macroscopic grids are very large, computers can be unstable by writing all information of all macroscopic grid points at the same time. To avoid that, the writings are divided by specifying this option. Writings will be done by each ``nmacro_write_group`` macroscopic grid points. (this number must be aliquot part of the total number of macroscopic grid points)


&maxwell
--------

- **al_em(3)** (real(8), Default=0d0)
   Size of simulation box in electromagnetic analysis. Unit of the length can be chosen by ``&units/unit_system``.

- **dl_em(3)** (real(8), Default=0d0)
   Spacing of real-space grids in electromagnetic analysis. Unit of length can be chosen by ``&units/unit_system``.

- **dt_em** (real(8), Default=0)
   Time step in electromagnetic analysis. Unit of time can be chosen by ``&units/unit_system``.

- **nt_em** (integer, Default=0)
   Number of total time steps for real-time propagation in electromagnetic analysis.

- **boundary_em(3,2)** (character, Default='default')
   Boundary condition in electromagnetic analysis. The first index(1-3 rows) corresponds to x, y, and z axes. The second index(1-2 columns) corresponds to bottom and top of the axes. If ``&system/yn_periodic='n'``, ``'default'``, ``'pml'``, and ``'pec'`` can be chosen. ``'pml'`` is absorbing boundary and ``'pec'`` is perfect electric conductor. ``'default'`` is ``'pml'``. If ``&system/yn_periodic='y'``, ``'default'``, ``'pml'``, and ``'periodic'`` can be chosen. ``'periodic'`` is periodic boundary. ``'default'`` is ``'periodic'``.

- **shape_file** (character, Default='none')
   Name of shape file in electromagnetic analysis. The shape files can be generated by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **media_num** (integer, Default=0)
   Number of media in electromagnetic analysis.

- **media_type(:)** (character, Default='vacuum')
   Type of media in electromagnetic analysis. ``'vacuum'``, ``'constant media'``, ``'pec'``, and ``'lorentz-drude'`` can be chosen. If ``'lorentz-drude'`` is chosen, linear response calculation can be done by ``&emfield/ae_shape1 or ae_shape2='impulse'``.

- **epsilon_em(:)** (real(8), Default=1d0)
   Relative permittivity of the media in electromagnetic analysis.

- **mu_em(:)** (real(8), Default=1d0)
   Relative permeability of the media in electromagnetic analysis.

- **sigma_em(:)** (real(8), Default=0d0)
   Conductivity of the media in electromagnetic analysis.

- **pole_num_ld(:)** (integer, Default=1)
   Number of poles of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis.

- **omega_p_ld(:)** (real(8), Default=0d0)
   Plasma frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis.
- **f_ld(:,:)** (real(8), Default=0d0)
   Oscillator strength of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``media_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``.

- **gamma_ld(:,:)** (real(8), Default=0d0)
   Collision frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``media_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``.

- **omega_ld(:,:)** (real(8), Default=0d0)
   Oscillator frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``media_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``.

- **wave_input** (character, Default='none')
   If ``'source'``, the incident pulse in electromagnetic analysis is generated by the incident current source.

- **ek_dir1(3)/ek_dir2(3)** (real(8), Default=0d0)
   Propagation direction of the first/second pulse.

- **source_loc1(3)/source_loc2(3)** (real(8), Default=0d0)
   Location of the incident current source of the first/second pulse. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/yn_periodic='n'`` while ranges from ``0`` to ``al_em`` for ``&system/yn_periodic='y'``.

- **obs_num_em** (integer, Default=0)
   Number of observation point in electromagnetic analysis. From the obtained results, figure and animation files can be generated by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **obs_samp_em** (integer, Default=1)
   Sampling time-step of the observation in electromagnetic analysis.

- **obs_loc_em(:,3)** (real(8), Default=0d0)
   Location of the observation point in electromagnetic analysis. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/yn_periodic='n'`` while ranges from ``0`` to ``al_em`` for ``&system/yn_periodic='y'``.

- **yn_obs_plane_em(:)** (character, Default='n')
   Enable(``'y'``)/disable(``'n'``). Output of the electrmagnetic fields on the planes (xy, yz, and xz planes) for each observation point. This option must be ``'y'`` for generating animation files by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **yn_wf_em** (character, Default='y')
   Enable(``'y'``)/disable(``'n'``). Applying a window function for linear response calculation when ``&calculation/theory=Maxwell``.

&analysis
---------

- **projection_option** (character, Default='no')
   | Old infomation: 3d
   Methods of projection.
   
  - ``'no'``: no projection.
  - ``'gs'``: projection to eigenstates of ground-state Hamiltonian.
  - ``'rt'``: projection to eigenstates of instantaneous Hamiltonian.
  

- **projection_decomp** (character, Default='n')[Trial]
   | Old infomation: 3d
   If ``'atom'`` combined with ``projection_option='gs'``, 
   the number of excited electron is decomposed into each atom 
   (this is printed in ``SYSname``\_nex_atom.data).

- **nenergy** (integer, Default=1000)
   Number of energy grids for frequency-domain analysis.
   This parameter is required when `'impulse'` is choosen in `&emfield/ae_shape1|2`.

- **de** (real(8), Default=0.01d0 eV)
   Energy spacing for analysis.
   Unit of energy can be chosen by ``&units/unit_energy``
   This parameter is required when `'impulse'` is choosen in `&emfield/ae_shape1|2`.

- **yn_out_psi** (character, Default='n')
   If ``'y'``, wavefunctions are output.
   For periodic system (``yn_periodic='y'``), it works only for ground state calculation. The converged wave functions of all orbitals with all k-points are printed in gs_wfn_cube or gs_wfn_vtk directory. The format is speficied by ``format3d``. 

- **yn_out_dos** (character, Default='n')
   If ``'y'``, density of state is output.

- **yn_out_dos_set_fe_origin** (character, Default='n')
   If ``'y'``, the electron energy is shifted to fix the Fermi energy as zero point.
   For ``&system/yn_periodic`` is ``'n'``, `` out_dos_fshift`` is not used 
   if ``&system/nstate`` is equal to ``&system/nelec``/2.

- **out_dos_start** (real(8), Default=-1d10 eV)
   Lower bound (energy) of the density of state spectra.
   If this value is lower than a specific value near the lowest energy level, 
   this value is overwritten by that value. 

- **out_dos_end** (real(8), Default=1d10 eV)
   Upper bound (energy) of the density of state spectra.
   If this value is higher than a specific value near the highest energy level, 
   this value is overwritten by that value. 

- **out_dos_nenergy** (integer, Default=601)
   Number of  energy points sampled in the density of state spectra.
 
- **out_dos_width** (real(8), Default=0.1d0 eV)
   Smearing width used in the density of state spectra..

- **out_dos_function** (character, Default='gaussian')
   Choise of smearing method for the density of state spectra..
   ``gaussian`` and ``lorentzian`` function are available.

- **yn_out_pdos** (character, Default='n')
   | Old infomation: 0d
   If ``'y'``, projected density of state is output.

- **yn_out_dns** (character, Default='n')
   If ``'y'``, the spatial electron density distribution at the ground state is output.

- **yn_out_dns_rt/out_dns_rt_step** (Character/Integer, Default='n')
   If ``'y'``,  the spatiotemporal electron density distribution during real-time time-propagation is output
   every ``outdns_rt_step`` time steps.

- **yn_out_dns_trans/out_dns_trans_energy** (Character/Real(8), Default='n'/1.55d0eV)[Trial]
   | Old infomation: 3d
   If ``'y'``, transition in different density from the ground state at specified field frequency omega(given by ``out_dns_trans_energy``) is calculated by drho(r,omega)=FT(rho(r,t)-rho_gs(r))/T.

- **yn_out_elf** (character, Default='n')
   | Old infomation: 0d
   If ``'y'``, electron localization function is output.

- **yn_out_elf_rt/out_elf_rt_step** (Character/Integer,Default='n'/50)
   | Old infomation: 0d
   If ``'y'``, electron localization function 
   during real-time time-propagation is output
   every ``out_elf_rt_step`` time steps.

- **yn_out_estatic_rt/out_estatic_rt_step** (Character/Integer, Default='n'/50)
   | Old infomation: 0d
   If ``'y'``, static electric field
   during real-time time-propagation is output
   every ``out_estatic_rt_step`` time steps.

- **yn_out_rvf_rt/out_rvf_rt_step** (Character/Integer, Default='n'/10)[Trial]
   | Old infomation: 3d
   If ``'y'``, coordinates[A], velocities[au], forces[au] on atoms
   during real-time time-propagation are printed in ``SYSname``\_trj.xyz
   every ``out_rvf_rt_step`` time steps.
   If ``use_ehrenfest_md='y'``, 
   the printing option is automatically turned on.
   
- **yn_out_tm** (character, Default='n')[Trial]
   | Old infomation: 3d
   If ``'y'``, transition moments between occupied and virtual orbitals are printed into ``SYSname``\_tm.data after the ground state calculation.

- **out_projection_step** (integer, Default=100)
   | Old infomation: 3d
   Interval time step of projection analysis 
   if ``projection_option`` is not ``'no'``.
   
- **out_ms_step** (integer, Default=100)
   xxx.

- **format_voxel_data** (character, Default='cube')
   File format for three-dimensional volumetric data.
   ``'avs'``, ``'cube'``, and ``'vtk'`` can be chosen.

- **nsplit_voxel_data** (integer, Default=1)
   | Old infomation: 0d
   Number of separated files for three dimensional data.
   Effective only when ``format3d`` is ``'avs'``.
   ``numfiles_out_3d`` must be less than or equal to number of processes.

- **timer_process** (character, Default='n')[Trial]
   | Old infomation: 0d
   Basically, elapsed times are written in the output file. 
   But if ``timer_process`` is ``'y'``, 
   files of elapsed times for every process are also generated. 
   This variable is effective only for the real-time caululation.


&poisson
--------

- **layout_multipole** (character, Default=3)
   | Old infomation: 0d
   A variable to determine how to put multipoles in the Hartree potential calculation.

  - ``1``: A single pole is put at the center.
  - ``2``: Multipoles are put at the center of atoms.
  - ``3``: Multipoles are put at the center of mass of electrons in prepared cuboids.

- **num_multipole_xyz(3)** (integer, Default=0)
   | Old infomation: 0d
   Number of multipoles when ``meo`` is ``3``. When default is set, number of multipoles is calculated automatically.

- **threshold_cg** (real(8), Default=1d-15 a.u.)
   xxx.


&ewald
------

- **newald** (integer, Default=4)
   | Old infomation: 3d
   Parameter for Ewald method. 
   Short-range part of Ewald sum is calculated within ``newald`` th
   nearlist neighbor cells.

- **aewald** (real(8), Default=0.5d0)
   | Old infomation: 3d
   Square of range separation parameter for Ewald method in atomic unit. 


&opt[Trial]
-------------

- **nopt** (integer, Default=100)
   xxx.

- **convrg_opt_fmax** (real(8), Default=1d-3)[Trial]
   | Old infomation: 3d
   Convergence threshold of optimization in maximum force.


&md[Trial]
-----------
- **ensemble** (character, Default='NVE')[Trial]
   | Old infomation: 3d
   Ensemble in MD option: "NVE" or "NVT".

- **thermostat** (character, Default='nose-hoover')[Trial]
   | Old infomation: 3d
   Thermostat in "NVT" option: (currently only ``nose-hoover``).

- **step_velocity_scaling** (integer, Default=-1)[Trial]
   | Old infomation: 3d
   Time step interval for velocity-scaling. Velocity-scaling is applied if this is set to positive.

- **step_update_ps/step_update_ps2** (Integer/Integer, Default=10/1)[Trial]
   | Old infomation: 3d
   Time step interval for updating pseudopotential (Larger number makes calculation time reduce greatly, but gets inaccurate) in case of ``use_ehrenfest_md=y``. ``step_update_ps`` is for full update and ``step_update_ps2`` is for update without changing grid points array.

- **temperature0_ion_k** (real(8), Default=298.15d0)[Trial]
   | Old infomation: 3d
   Setting temperature [K] for NVT ensemble, velocity scaling and generating initial velocities.

- **yn_set_ini_velocity** (character, Default='n')[Trial]
   | Old infomation: 3d
   Initial velocities are set.

  - ``y``: Generate initial velocity with Maxwell-Bortzman distribution.
  - ``r``: Read initial velocity from file specified by keyword of ``file_ini_velocity``. This is, for example, used for restarting MD from the previous run. The last atomic coordinates and velocities are printed in ``SYSname``\_trj.xyz. (atomic coordinate also should be copied from the previous output and put in the next input file for restart)
    
- **file_ini_velocity** (character, Default='none')[Trial]
   | Old infomation: 3d
   File name for initial velocities. This is read when ``set_ini_velocity`` is ``'r'``. The format is simply vx(iatom) vy(iatom) vz(iatom) in each line. The order of atoms must be the same as the given coordinates in the main input file. In case of using nose-hoover thermostat, a thermostat variable should be put at the last line (all atomic unit). 

- **thermostat_tau** (real(8), Default=41.34d0 a.u. or 1d0 fs)[Trial]
   | Old infomation: 3d
   Parameter in Nose-Hoover method: controlling time constant for temperature.
   Default is ``41.34[au] or 1.0[fs]``.

- **seed_ini_velocity** (integer, Default=123)[Trial]
   | Old infomation: 3d
   Random seed (integer number) to generate initial velocity if ``set_ini_velocity`` is set to y.
   Default is ``123``.

- **yn_stop_system_momt** (character, Default='n')[Trial]
   | Old infomation: 3d
   Center of mass is stopped every time step.


&code
-----

- **yn_want_stencil_openmp_parallelization** (character, Default='y')[Trial]
   xxx.

- **yn_want_stencil_hand_vectorization** (character, Default='y')[Trial]
   xxx.

- **yn_force_stencil_openmp_parallelization** (character, Default='n')[Trial]
   xxx.

- **yn_force_stencil_sequential_computation** (character, Default='n')[Trial]
   xxx.

- **yn_want_communication_overlapping** (character, Default='n')[Trial]
   xxx.
   

**Following variables are moved from the isolated part. Some of them may be added to common input, be combined to it, and be removed.**


&group_fundamental[Trial]
-------------------------

- **iditer_nosubspace_diag** (integer, Default=10)[Trial]
   | Old infomation: 0d
   Iterations for which subspace diagonalization is not done if ``&scf/subspace_diagonalization`` is ``'y'``.
   
- **ntmg** (integer, Default=1)[Trial]
   | Old infomation: 0d
   Number of multigrid calculation for gs. At the moment, there is a malfunction in this variable, and recovery is needed.

- **idisnum(2)** (integer, Default=1,2)[Trial]
   | Old infomation: 0d
   Label numbers for two atoms which are measured the distance. 

- **iwrite_projection** (integer, Default=0)[Trial]
   | Old infomation: 0d
   A variable for projection. 

- **itwproj** (integer, Default=-1)[Trial]
   | Old infomation: 0d
   The projection is calculated every ``itwproj`` time steps. 

- **iwrite_projnum** (integer, Default=0)[Trial]
   | Old infomation: 0d
   There is a malfunction in this variable.

- **itcalc_ene** (integer, Default=10)[Trial]
   | Old infomation: 0d
   Total energy is calculated every ``itcalc_ene`` time steps. There may be a malfunction in this variable.


&group_hartree[Trial]
----------------------

- **lmax_lmp** (integer, Default=4)[Trial]
   | Old infomation: 0d
   A maximum angular momentum for multipole expansion in the Hartree-cg calculation. 


&group_others[Trial]
---------------------

- **iswitch_orbital_mesh** (integer, Default=0)[Trial]
   | Old infomation: 0d
   A variable to apply descending order for orbitals in the ground state calculation.

- **iflag_psicube** (integer, Default=0)[Trial]
   | Old infomation: 0d
   A variable to generate cube files for wave functions. This variable will be removed.

- **num_projection** (Interger, Default=1)[Trial]
   | Old infomation: 0d
   Number of orbitals for projections.

- **iwrite_projection_ob(200)** (Interger, Default=1, 2, 3, ..., 200)[Trial]
   | Old infomation: 0d
   Orbital number to be written as projections.

- **iwrite_projection_k(200)** (Interger, Default=1)[Trial]
   | Old infomation: 0d
   This variable will be removed.

- **filename_pot** (character, Default='pot')[Trial]
   | Old infomation: 0d
   Name of file to be written local potentials. 

- **iwrite_external** (integer, Default=0)[Trial]
   | Old infomation: 0d
   A variable to generate file to be written local potentials. 

- **iflag_intelectron** (integer, Default=0)[Trial]
   | Old infomation: 0d
   A variable related to the quadrupole caluclation.

- **num_dip2** (integer, Default=1)[Trial]
   | Old infomation: 0d
   Number of area where dipole moments are calculated.

- **dip2boundary(100)** (real(8), Default=0d0 a.u.)[Trial]
   | Old infomation: 0d
   Boundary position of area where dipole moments are calculated.

- **dip2center(100)** (real(8), Default=0d0 a.u.)[Trial]
   | Old infomation: 0d
   Origin in the dipole moment calculation. 

- **itotntime2** (integer, Default=0)[Trial]
   | Old infomation: 0d
   Number of time steps in the reentrance for real-time calculation.
   There may be a malfunction in this variable.

- **iwdenoption** (integer, Default=0)[Trial]
   | Old infomation: 0d
   A variable to determine whether 3d output is generated in real-time calculation. 
   This variable will be removed.

- **iwdenstep** (integer, Default=0)[Trial]
   | Old infomation: 0d
   3d output is generated every ``iwdenstep`` time steps.
   This variable will be removed.

- **iflag_estatic** (integer, Default=0)[Trial]
   | Old infomation: 0d
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



