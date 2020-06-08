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
-  `&group_hartree[Trial]`_ 
-  `&group_others[Trial]`_


&calculation
------------

- **theory** (character, default='none')
   | Choice of Calculation theories.
   | Options
   |   ``dft``  / ground state calculation based on DFT
   |   ``dft_md``  / adiabatic ab initio MD simulations based on DFT
   |   ``tddft_response``  / simulations based on TDDFT for response
   |   ``tddft_pulse``  / simulations based on TDDFT using pulsed light
   |   ``single_scale_maxwell_tddft``  / coupled Maxwell and TDDFT single-scale simulation
   |   ``multi_scale_maxwell_tddft``  / coupled Maxwell and TDDFT multi-scale simulation
   |   ``maxwell``  / electromagnetic simulations based on the Maxwell's equations
   |   ``dft_k_expand`` / convert checkpoint data of dft with k-points calculation to that of larger supercell system with gamma-point

- **yn_md** (character, Default='n')[Trial]
   | Available for ``theory='dft'`` (Adiabatic ground-state MD) and ``theory='tddft_pulse'`` (Ehrenfest MD).
   | Molecular dynamics option.
   Enable(``'y'``)/disable(``'n'``). 

- **yn_opt** (character, Default='n')[Trial]
   | Available for ``theory='dft'``.
   | Geometry optimization option.
   Enable(``'y'``)/disable(``'n'``).

&control
--------

- **sysname** (character, Default='default')
   | Available for all options of ``theory``.
   Name of calculation. This is used for a prefix of output files.

- **base_directory** (character, Default='./')
   | Available for all options of ``theory``.
   Name of a default directory, where the basic results will be written down.

- **yn_restart** (character, Default='n')
   | Available for ``theory='dft*' or '*tddft*'``.
   | Restart option.
   Enable(``'y'``)/disable(``'n'``).

- **directory_read_data** (character, Default='restart/')
   | Directory name for the restart data that is written down in the previous run 

- **yn_self_checkpoint** (character, Default='n')
   | If set `'y'`: When saving intermediate results of the simulation (this call checkpointing), each process write/read a checkpoint data independently.
   | This option helps large-scale simulation to recover from system failure, which reduce restart costs.

- **checkpoint_interval** (integer, Default=0)
   | Available for ``theory='dft*' or '*tddft*'``.
   Interval of time step (or iteration step) of writing down check-point data during the time-propagation or iteration. 
   These are not written down If ``0`` is set.

- **yn_reset_step_restart** (character, Default='n')
   | Available for ``yn_restart='y'`` with the DFT/TDDFT based options of ``theory``.
   In the case of restarting, the initial step of SCF iteration (for DFT) or time step (for TDDFT) are reset to 0 at begining.
   Then, the memory of the density in the previous SCF iteration steps (in GS) is abondoned.

- **read_gs_restart_data** (character, Default='all')
   | Available for ``yn_restart='y'`` with ``theory='dft'``.
   | Options
   |   ``all``  / all of restart data are read
   |   ``rho_inout``  / only electron densities including memories at previous iteration steps are read
   |   ``rho``  / only the latest electron density is read (user-made data)
   |   ``wfn``  / only wavefunctions is read
   Specified data which is included in the restart (or checkpoint) directory generated in the previous calculation is used for restarting SCF iteration in DFT.
   The default option ``'all'`` gives the complete restart. The other options use a part of restart data (other necessary data is generated as done in the initial SCF step)

- **write_gs_restart_data** (character, Default='all')
   | Available for ``theory='dft'``.
   | Options
   |   ``all``  / all of restart data are written out
   |   ``rho_inout``  / only electron densities including memories at previous iteration steps are written out
   |   ``wfn``  / only wavefunctions is written out
   Specified data is written out in the restart (or checkpoint) directory. 
   The default option ``'all'`` gives the complete set of restart data. 


- **time_shutdown** (real(8), Default=-1d0)[Trial]
   | Available for ``theory='XXX'``.
   Timer for automatic shutdown. The unit is second.
   If negative time is chosen, the automatic shutdown is not performed.

- **method_wf_distributor** (character, Default='single')
   | Available for ``theory='dft*' or '*tddft*'``.
   | Select a method of save/load the wave function.
   | 'single': wave function saves/loads a single shared file.
   | 'slice' : wave function saves/loads to a file per the orbital function. 
   | 'slice' reduces I/O costs, and they can helps flexible large-scale simulation.

- **nblock_wf_distribute** (integer, Default='16')
   | Available for ``method_wf_distributor='slice'``.
   | 'slice' mode saves ``nblock_wf_distribute``-files to a directory.
   | In a default, they will saves 16 files to same directory.

- **yn_gbp** (character, Default='n')
   | XXX.

- **yn_gbp_fourier0** (character, Default='n')
   | XXX

&units
------

- **unit_system** (character, Default='au')
   | Units of input variables. 
   | Options
   |   ``'au'`` or ``'a.u.'`` / atomic unit system.
   |   ``'A_eV_fs'`` / Angstrom-eV-fs unit system

&parallel
---------

- **nproc_k/nproc_ob/nproc_rgrid(3)** (integer, Default=0)
   | Old infomation: 0d
   | Options
   |   ``nproc_k``/ Number of MPI parallelization for orbitals that related to the wavefunction calculation.
   |   ``nproc_ob``/ Number of MPI parallelization for orbitals that related to the wavefunction calculation.
   |   ``nproc_rgrid(3)'``/ Number of MPI parallelization for each direction in real-space that related to the wavefunction and the electron density calculations. 
   |
   | Defaults are ``0`` for ``nproc_k``/``nproc_ob`` and ``(0,0,0)`` for ``nproc_rgrid``. If users use the defaults, automatic proccess assignment is done. Users can also specify ``nproc_k``, ``nproc_ob``, and ``nproc_rgrid`` manually. In that case, ``nproc_k`` must be set to ``1`` for isolated system calculations. ``nproc_k`` and ``nproc_k`` must be set to ``1`` for ``theory='maxwell'``. In addition, followings must be satisfied.
   |
   |   ``nproc_k`` \* ``nproc_ob`` \* ``nproc_rgrid(1)`` \* ``nproc_rgrid(2)`` \* ``nproc_rgrid(3)`` \= total number of processes.

- **yn_ffte** (character, Default='n')
   | Available for ``&system/yn_periodic='y'``
   | Old infomation: 0d
   | Method of Fourier transformation.  
   | Enable(``'y'``)/disable(``'n'``).
   | SALMON uses FFT (via FFTE library) to solve poisson equation.
   | When enabling it, followings must be satisfied.
   |
   | ``mod(num_rgrid(1), nproc_rgrid(2)) == 0``
   | ``mod(num_rgrid(2), nproc_rgrid(2)) == 0``
   | ``mod(num_rgrid(2), nproc_rgrid(3)) == 0``
   | ``mod(num_rgrid(3), nproc_rgrid(3)) == 0``

- **yn_scalapack** (character, Default='n')
   | Available for ``&calculation/theory='dft' or 'dft_md'``
   | SALMON uses ScaLAPACK library to solve eigenvalue problem in subspace diagonalization.
   | When enabling it, you should build SALMON by linking ScaLAPACK library.

- **yn_eigenexa** (character, Default='n')
   | Available for ``&calculation/theory='dft' or 'dft_md'``
   | SALMON uses RIKEN R-CCS EigenExa library to solve eigenvalue problem in subspace diagonalization.
   | When enabling it, you should build SALMON by linking ScaLAPACK and EigenExa libraries.
   
- **yn_scalapack_red_mem** (character, Default='n')
   | Available for ``&parallel/yn_scalapack='y'`` or ``&parallel/yn_eigenexa='y'``
   | We use ScaLAPACK/EigenExa libraries by optimized algorithm to reduce memory consumption.

- **process_allocation** (character, Default='grid_sequential')
   | Old infomation: 0d
   | You can select the process allocation ordering.
   | ``'grid_sequential'``    / real-space grid major ordering.
   | ``'orbital_sequential'`` / orbital-space major ordering.
   |
   | Suggestion:
   |   ``&calculation/theory='dft' or 'dft_md'``            / `orbital_sequential`
   |   ``&calculation/theory='tddft*' or '*maxwell_tddft'`` / `grid_sequential`


&system 
-------

- **yn_periodic** (character, Default='n')
   | Available for all options of ``theory``.
   | Option of periodic boundary condition.
   |   ``'y'`` / periodic systems (solids)
   |   ``'n'`` / isolated systems

- **spin** (character, Default='unpolarized')
   | Available for all options of ``theory`` except for ``theory='maxwell'``.
   | Variable for classification of spin-unpolarized (closed shell) systems and spin-polarized (open shell) systems. 
   | Options 
   |   ``'unpolarized'`` / spin-unpolarized systems (default)
   |   ``'polarized'`` / spin-polarized systems

- **al(3)** (real(8), Default=0d0)
   | Available for all options of ``theory`` except for ``theory='maxwell'``.
   Spatial grid box size or lattice constants for cuboid cell (x, y, z).
   For nonorthogonal cell, see ``al_vec1,al_vec2,al_vec3``.

- **al_vec1(3)/al_vec2(3)/al_vec3(3)** (real(8), Default=0d0)
   | Available for all options of ``theory`` except for ``theory='maxwell'``.
   Primitive lattice vectors for nonorthogonal cell.

- **isym** (integer, Default=1)
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Number of symmetries that can be used for reduction of k-points.

- **crystal_structure** (character, Default='none')
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Name of symmetry that can be used for the reduction of # of k-points.

- **nstate** (integer, Default=0)
   | Available for the DFT/TDDFT based options of ``theory``.
   Number of orbitals/bands.

- **nstate_spin(2)** (integer, Default=0)
   | Available for the DFT/TDDFT based options of ``theory``.
   Number of orbitals/bands for up/down-spin electrons can be specified for each by ``nstate_spin(1)/nstate_spin(2)``.
   This option is incompatible with ``nstate``  (?? does it mean nstate specified is ignored if this option is specified ??)

- **nelec** (integer, Default=0)
   | Available for the DFT/TDDFT based options of ``theory``.
   Number of valence electrons.

- **nelec_spin(2)** (integer, Default=0)
   | Available for the DFT/TDDFT based options of ``theory``.
   Number of up/down-spin electrons can be specified for each by ``nelec_spin(1)/nelec_spin(2)``.
   This option is incompatible with ``nelec`` (?? does it mean nelec specified is ignored if this option is specified ??)

- **temperature** (real(8), Default=-1d0)
   | Available for DFT-based options of ``theory`` 
   Temperature of electrons. When you calculate a system of zero band-gap energy like metals, zero or positive number of the temperature should be given (Default is ``-1.0`` : this is for system which has a band gap energy). The value must be given by the unit of energy as specified in ``&units/unit_system``. The kelvin unit can be used by the keyword ``temperature_k`` (see next). 
   
- **temperature_k** (real(8), Default=-1d0)[Trial]
   | Available for DFT-based options of ``theory``
   The same as ``temperature`` but in kelvin.

- **nelem** (integer, Default=0)
   | Available for the DFT/TDDFT based options of ``theory``.
   Number of used atomic elements in the system.

- **natom** (integer, Default=0)
   | Available for the DFT/TDDFT based options of ``theory``.
   Number of atoms in the system.

- **file_atom_red_coor** (character, Default='none')[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   File name for atomic positions given in reduced coordinates. 
   This option is incompatible with ``&system/file_atom_coor``, ``&atomic_coor``, and ``&atomic_red_coor``.

- **file_atom_coor** (character, Default='none')[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   File name for atomic Cartesian coordinates (The unit is specified by ``&units/unit_system``). 
   This option is incompatible with ``&system/file_atom_coor``, ``&atomic_coor``, and ``&atomic_red_coor``.
   (XXX why this keyword is not in &atomic_coor ?? XXX)

&atomic_red_coor
----------------

Positions of atoms is given in reduced coordinates as follows:

|  'Si'	 0.00  0.00  0.00  1
|  'Si'	 0.25  0.25  0.25  1
|  ...

Here, the information of atoms is ordered in row. 
For example, the first row is for the first atom. 
The number of rows must be equal to ``&system/natom``.
The first coloum can be any caracters and does not affect calculations.
The second, third and fourth columns are reduced coordinates for
the first, second and third directions, respectively. 
The fifth column is a serial number of the atom spieces, which is defined in ``&pseudo``.
This option is incompatible with 
``&system/file_atom_red_coor``, ``&system/file_atom_coor``, and ``&atomic_coor``.


&atomic_coor
------------

Cartesian atomic coordinates.
The format is the same as &atomic_red_coor.
The unit can be chosen by ``&units/unit_length``.
This option is incompatible with 
``&system/file_atom_red_coor``, ``&system/file_atom_coor``, and ``&atomic_red_coor``.


&pseudo
-------

Input for psudopotentials. Size of array (:) is equal to ``&system/nelem``.

- **izatom(:)** (integer, Default=-1)
   | Available for the DFT/TDDFT based options of ``theory``.
   Atomic number.

- **file_pseudo(:)** (character, Default='none')
   | Available for the DFT/TDDFT based options of ``theory``.
   File name for pseudopotential.

- **lmax_ps(:)** (integer, Default=-1)
   | Available for the DFT/TDDFT based options of ``theory``.
   Maximum angular momentum of pseudopotential projectors. If not given, it is automatically read from the pseudopotential file.

- **lloc_ps(:)** (integer, Default=-1)
   | Available for the DFT/TDDFT based options of ``theory``.
   Angular momentum of pseudopotential that will be treated as local.

- **yn_psmask(:)** (character, Default='n')[Trial]
   | Available for the DFT/TDDFT based options of ``theory``.
   | Fourier filtering for pseudopotentials. 
   Enable(``'y'``)/disable(``'n'``) 

- **alpha_mask(:)** (real(8), Default=0.8d0)[Trial]
   | Available for the DFT/TDDFT based options of ``theory``.
   Parameter for the Fourier filtering for pseudopotential.

- **gamma_mask(:)** (real(8), Default=1.8d0)[Trial]
   | Available for the DFT/TDDFT based options of ``theory``.
   Parameter for the Fourier filtering for pseudopotential.

- **eta_mask(:)** (real(8), Default=15.0d0)[Trial]
   | Available for the DFT/TDDFT based options of ``theory``.
   Parameter for the Fourier filtering for pseudopotential.


&functional
-----------

- **xc** (character, Default='none')
   | Available for the DFT/TDDFT based options of ``theory``.
   | Exchange-correlation functionals.
   | At present version, the functional 'PZ', 'PZM' and 'TBmBJ' is available for both 0d/3d calculations, and the functionals 'TPSS' and 'VS98' are available for 3d calculations. (XXX need check XXX)
   | Options
   |   ``'PZ'``: Perdew-Zunger LDA :Phys. Rev. B 23, 5048 (1981).
   |   ``'PZM'``: Perdew-Zunger LDA with modification to improve sooth connection between high density form and low density one. :J. P. Perdew and Alex Zunger, Phys. Rev. B 23, 5048 (1981).
   |   ``'TBmBJ'``: Tran-Blaha meta-GGA exchange with Perdew-Wang correlation. :Fabien Tran and Peter Blaha, Phys. Rev. Lett. 102, 226401 (2008). John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992).
   |   ``'TPSS'``: Tao, Perdew, Staroverov and Scuseria meta-GGA exchange correlation. :J. Tao, J. P. Perdew, V. N. Staroverov, and G. E. Scuseria, Phys. Rev. Lett. 91, 146401 (2003).
   |   ``'VS98'``:  van Voorhis and Scuseria exchange with Perdew-Wang correlation: T. Van Voorhis and G. E. Scuseria, J. Chem. Phys. 109, 400 (1998).

- **cname, xname** (character, Default='none')
   | Available for ``theory='XXX'``.
   | XXX

- **alibxc, alibx, alibc** (character, Default='none')
   | Available for the DFT/TDDFT based options of ``theory``.
   By specifying ``alibxc``, the functionals prepared in libxc package are available. 
   They can be set indivisually by specifying ``alibx`` and ``alibc``.
   To use libxc libraries, ``--with-libxc`` option must be added in excecuting configure. 
   The available option of the exchange-correlation functionals are listed in the LibXC website. 
   [See http://www.tddft.org/programs/libxc/functionals/]
   
- **cval** (real(8), Default=-1d0)
   | Available for ``xc='TBmBJ'``.
   Mixing parameter in Tran-Blaha meta-GGA exchange potential. If ``cval`` is set to a minus value, the mixing-parameter computed
   by the formula in the original paper [Phys. Rev. Lett. 102, 226401 (2008)].
   Default is estimated from :math:`\left\langle |\nabla \rho(\mathbf{r};t)| / \rho(\mathbf{r};t) \right\rangle`.


&rgrid
------

- **dl(3)** (real(8), Default=0d0)
   | Available for the DFT/TDDFT based options of ``theory``.
   Spacing of real-space grids. 
   (This cannot be used together with ``&rgrid/num_rgrid``.)
   If ``&system/yn_periodic='y'``,
   the grid spacing is automatically adjusted in calculations
   so that the grid box size ``&system/al(3)`` becomes divisible by the spacing.  (XXX really?? XXX)

- **num_rgrid(3)** (integer, Default=0)
   | Available for the DFT/TDDFT based options of ``theory``.
   | Old infomation: 3d
   Dividing number of real-space grids for each direction.
   (This cannot be used together with ``&rgrid/dl``.)

&kgrid
------

- **num_kgrid(3)** (integer, Default=1)
   | Available for ``yn_periodic='y'``.
   Number of k-points (grid points of k-vector) for each direction discretizing the Brillouin zone.

- **file_kw** (character, Default='none')
   | Available for ``yn_periodic='y'``.
   File name for user specified k-points.
   This file will be read if ``num_kgrid`` is smaller than 1.


&tgrid
------

- **nt** (integer, Default=0)
   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   Number of total time steps for real-time propagation.

- **dt** (real(8), Default=0d0)
   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   Time step size.

- **gram_schmidt_interval** (integer, Default=0)
   | Available for TDDFT-based options of ``theory``.
   Interval of time step for the Gram-Schmidt orthonormalization of the orbital wavefunctions in the time-evolution calculation.

&propagation
------------

- **n_hamil** (integer, Default=4)[Trial]
   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   | Old infomation: 0d
   Order of Taylor expansion of a propagation operator.

- **propagator** (character, Default=middlepoint')
   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   | Propagator (time-integrator).
   | Options
   |   ``middlepoint`` / propagator with the Hamiltoinan at midpoint of two-times.
   |   ``etrs`` / time-reversal symmetry propagator.
   [M.A.L. Marques, A. Castro, G.F. Bertsch, and A. Rubio, Comput. Phys. Commun., 151 60 (2003)].

- **yn_fix_func** (character(1), Default='n')[Trial]
   | Available for 'dft_md' and TDDFT-based options of ``theory``.
   Option not to update functional (or Hamiltonian) in time-evolution, i.e., keep ground state Hamiltonian. (currently not available)

&scf
----

- **method_init_wf** (character, Default='gauss')
   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Options
   |   ``gauss`` / put single gauss function using a random number on each initial orbital
   |   ``gauss2`` / put two gauss functions using a random number on each initial orbital
   |   ``gauss3`` / put three gauss functions using a random number on each initial orbital
   |   ``gauss4`` / put four gauss functions using a random number on each initial orbital
   |   ``gauss5`` / put five gauss functions using a random number on each initial orbital
   |   ``gauss10`` / put ten gauss functions using a random number on each initial orbital
   |   ``random`` / give a random number at each real-space grid point on each initial orbital
   The generation method of the initial wavefunction (orbital) at the begening of the SCF iteration in DFT calculation.


- **iseed_number_change** (integer, Default=0)
   | Available for 'dft' and 'dft_md' options of ``theory``.
   The seed of the random numbers are changed by adding the given number for generating the initial wavefunctions.

- **nscf** (integer, Default=300)
   | Available for 'dft' and 'dft_md' options of ``theory``.
   Number of maximum SCF cycle in DFT calculation. 

- **method_min** (character, Default='cg') 
   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Method for SCF iteration
   | Options
   |  ``cg`` / Conjugate-Gradient(CG) method
   |  ``diis`` / DIIS method
   |  ``cg-diis`` / CG-DIIS method 

- **ncg** (integer, Default=4)
   | Available for 'dft' and 'dft_md' options of ``theory``.
   Number of interation of Conjugate-Gradient method for each scf-cycle.

- **ncg_init** (integer, Default=4)
   | Available for 'dft' and 'dft_md' options of ``theory``.
   Number of interation of Conjugate-Gradient method for the first SCF step.

- **method_mixing** (character, Default='broyden') 
   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Methods for density/potential mixing for scf cycle. ``simple`` and ``broyden`` can be chosen.
   | Options
   |  ``simple`` / Simple mixing method
   |  ``broyden`` / modified-Broyden method

- **mixrate** (real(8), Default=0.5d0)
   | Available for ``method_mixing='simple'`` in 'dft' and 'dft_md' options of ``theory``.
   | Mixing ratio for simple mixing.

- **nmemory_mb** (integer, Default=8)
   | Available for ``method_mixing='broyden'`` in 'dft' and 'dft_md' options of ``theory``.
   Number of previous densities to be stored in SCF iteration cycle for the modified-Broyden method. 
   If ``&system/yn_periodic`` is ``'n'``, ``nmemory_mb`` must be less than 21.

- **alpha_mb** (real(8), Default=0.75d0)
   | Available for ``method_mixing='broyden'`` in 'dft' and 'dft_md' options of ``theory``.
   Parameter of the modified-Broyden method.

- **nmemory_p** (integer, Default=4)
   | XXX

- **beta_p** (real(8), Default=0.75d0)
   | XXX

- **yn_auto_mixing** (character, Default='n')
   | XXX

- **update_mixing_ratio** (real(8), Default=3.0d0)
   | XXX

- **yn_subspace_diagonalization** (character, Default='y')
   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Option of subspace diagonalization during SCF cycle.
   | Options
   |   ``'y'`` / enable
   |   ``'n'`` / disable

- **convergence** (character, Default='rho_dne')
   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Quantity that is used for convergence check in SCF calculation. 
   | Options
   |   ``'rho_dne'``/ Convergence is checked by sum_ix|rho(ix,iter)-rho(ix,iter-1)|dx/N, where iter is iteration number of SCF calculation and N is ``&system/nelec``, the number of the valence electrons.
   |   ``'norm_rho'``/ Convergence is checked by the square of the norm of difference of density, ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\=sum_ix|rho(ix,iter)-rho(ix,iter-1)|\ :sup:`2`\. 
   |   ``'norm_rho_dng'``/ Convergence is checked by ||rho_iter(ix)-rho_iter-1(ix)||\ :sup:`2`\/(number of grids). "dng" means "devided by number of grids".
   |   ``'norm_pot'``/ Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\, where Vlocal is Vh + Vxc + Vps_local.
   |   ``'pot_dng'``/ Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||\ :sup:`2`\/(number of grids).

- **threshold** (real(8), Default=1d-17 [a.u.] (for ``convergence='rho_dne'``) and -1 (for other options of ``convergence``)
   | Available for 'dft' and 'dft_md' options of ``theory``.
   | Threshold for convergence that is specified by ``convergence`` keyword.
   | Unit conversions are: 1 a.u.= 45.54 A\ :sup:`-6`\ for ``convergence='norm_rho'`` and ``'norm_rho_dng'``,  1 a.u.= 33.72x10\ :sup:`4`\ A\ :sup:`-6`\eV\ :sup:`2`\ for  ``convergence='norm_pot'`` and ``'norm_pot_dng'``

- **iditer_notemperature** (integer, Default=10) 
   | Available for 'dft' and 'dft_md' options of ``theory``.
   | XXX

- **nscf_init_diagonal** (integer, Default=10)
   | Available for 'dft' option of ``theory``.
   The number of initial iterations for which subspace diagonalization is not done if ``&scf/yn_subspace_diagonalization`` is ``'y'``.

- **nscf_init_mix_zero** (Integer, Default=-1)
   | Available for 'dft' option of ``theory``.
   The densities is not mixed (i.e. fixed) during the given number of the SCF iteration cycle, that is, wavefunctions are optimized without updating the density. 

- **conv_gap_mix_zero** (real(8), Default=99999d0)
   | Available for positive number of ``nscf_init_mix_zero`` with 'dft' option of ``theory``.
   The condition to quite the fixed density forced by ``step_initial_mix_zero`` option.
   The density is allowed to start mixing after the band-gap energy exceeds the given number for consecutive five SCF iteration steps, 


&emfield
--------

- **trans_longi** (character, Default='tr')
   | Available for ``yn_periodic='y'`` with 'maxwell' and TDDFT based options of ``theory``.
   | Boundary condition for fields on macro-scale in solid-state calculations.
   | Options
   |   ``'tr'`` / Transverse  
   |   ``'lo'`` / longitudinal

- **ae_shape1/ae_shape2** (character, Default='none')
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Envelope shape of the first/second pulse.
   | Options
   |   ``'impulse'`` / Impulsive fields.
   |   ``'Acos2'`` / Envelope of cos\ :sup:`2`\ for a vector potential.
   |   ``'Acos3'`` / Envelope of cos\ :sup:`3`\ for a vector potential.
   |   ``'Acos4'`` / Envelope of cos\ :sup:`4`\ for a vector potential.
   |   ``'Acos6'`` / Envelope of cos\ :sup:`6`\ for a vector potential.
   |   ``'Acos8'`` / Envelope of cos\ :sup:`8`\ for a vector potential.   
   |   ``'Ecos2'`` / Envelope of cos\ :sup:`2`\ for a electric field.
   |   ``'Asin2cos'`` [Trial] / Envelope of sin\ :sup:`2`\cos for a vector potential.
   |   ``'Asin2cw'`` [Trial] /  XXX
   |   ``'Esin2sin'`` [Trial] / Envelope of sin\ :sup:`2`\sin for a electric field.
   |   ``'input'`` [Trial] / read-in user-defined field is used given by file name of 'input_Ac.dat'
   |   ``'none'`` can be also chosen.
   | XXX should be checked XXX

- **file_input1** (character, Default='')
   | XXX

- **e_impulse** (real(8), Default=1d-2 a.u.)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Momentum of impulsive perturbation. This valiable has the dimention of momentum, energy*time/length.

..
(commented out: not implemented yet)
- **t_impulse**
   | Available for ``theory='XXX'``.
   not yet implemented XXX
..
   
- **E_amplitude1/E_amplitude2** (real(8), Default=0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Maximum amplitude of electric fields for the first/second pulse.
   This valiable has the dimension of electric field, energy/(length*charge).
   This cannot be set with ``&emfield/I_wcm2_1`` (``I_wcm2_2``) simultaneously.

- **I_wcm2_1/I_wcm2_2** (real(8), Default=-1d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Peak intensity (W/cm\ :sup:`2`\) of the first/second pulse.
   This valiable cannot be set with ``&emfield/E_amplitude1`` (``E_amplitude2``) simultaneously.

- **tw1/tw2** (real(8), Default=0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Duration of the first/second pulse (edge-to-edge time length). 

- **omega1/omega2** (real(8), Default=0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   | Mean photon energy (average frequency multiplied by the Planck constant) of the first/second pulse. 

- **epdir_re1(3)/epdir_re2(3)** (real(8), Default=1d0, 0d0, 0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   Real part of polarization unit vector for the first/second pulse.

- **epdir_im1(3)/epdir_im2(3)** (real(8), Default=0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   Imaginary part of polarization unit vector for the first/second pulse.

- **phi_cep1/phi_cep2** (real(8), Default=0d0/0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   Carrier emvelope phase of the first/second pulse.

- **t1_t2** (real(8), Default=0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   Time-delay between the first and the second pulses.

- **t1_start** (real(8), Default=0d0)
   | Available for 'maxwell' and TDDFT based options of ``theory``.
   Time-delay of the first pulse. (this is not available for multiscale option).

- **num_dipole_source** (integer, Default=0)
   | Available for TDDFT based options of ``theory``.
   Number of radiation sources for optical near fields.
   Maximum number is ``2``.

- **vec_dipole_source(3,num_dipole_source)** (real(8), Default=0d0)
   | Available for TDDFT based options of ``theory``.
   Dipole vectors of the radiation sources for the optical near fields.

- **cood_dipole_source(3,num_dipole_source)** (real(8), Default=0d0)
   | Available for TDDFT based options of ``theory``.
   Central coordinates of the dipole vectors of the radiation sources.

- **rad_dipole_diele** (real(8), Default=2d0 a.u.)
   | Available for TDDFT based options of ``theory``.
   Radii of dielectric spheres for the radiation sources.

- **cutoff_G2_emfield** (real(8), Default=-1d0)
   | XXX


&multiscale
-----------

- **fdtddim** (character, Default='1d')[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Dimension of FDTD calculation for multi-scale Maxwell-Kohn-Sham method.

- **twod_shape** (character, Default='periodic')[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Boundary condision of the second dimension for FDTD calculation with 
   multi-scale Maxwell-Kohn-Sham method.

- **nx_m** (integer, Default=1)
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Number of macroscopic grid points inside materials for x-direction.

- **ny_m/nz_m** (integer, Default=1)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Number of macroscopic grid points inside materials for (y/z)-direction.

- **hx_m** (real(8), Default=0d0)
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Spacing of macroscopic grid points inside materials for (x)-direction.
   Unit of length can be chosen by ``&units/unit_length``.

- **hy_m/hz_m** (real(8), Default=0d0)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Spacing of macroscopic grid points inside materials for (y/z)-direction.
   Unit of length can be chosen by ``&units/unit_length``.

- **nxvacl_m/nxvacr_m** (integer, Default=1/0)
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Number of macroscopic grid points for vacumm region.
   ``nxvacl_m`` gives the number for negative x-direction in front of material,
   while ``nxvacr_m`` gives the number for positive x-direction behind the material.

- **nx_origin_m/ny_origin_m/nz_origin_m** (integer, Default=1)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Origin coordinat of the grid points.

- **file_macropoint** (character, Default='')[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   If file name is specified in the option, the coordinates of the macropoints are set from the file.

- **set_ini_coor_vel** (character, Default='n')[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   Set initial atomic coordinates and velocities for each macro-grid point. This must be given with specific directories and files: 
   Prepare ``directory``/multiscale/MXXXXXX/ini_coor_vel.dat, where 'XXXXXX' is the index number of the macro-grid point of the material region usually starting from '000001' up to the number of macro-grid point. The format of the file 'ini_coor_vel.dat' is just Rx, Ry, Rz, Vx, Vy, Vz (with space separation) for each atom (i.e. for each line), where the unit of the coordinates, Rx, Ry, Rz, is angstrom or a.u. speficied by ``unit_system`` but that of velocities is always a.u.. This option should be used together with ``read_gs_wfn_k_ms`` which is the option to read the ground state wave function for each macro-grid point. 

- **nmacro_write_group** (integer, Default=-1)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   If the number of macroscopic grids are very large, computers can be unstable by writing all information of all macroscopic grid points at the same time. To avoid that, the writings are divided by specifying this option. Writings will be done by each ``nmacro_write_group`` macroscopic grid points. (this number must be aliquot part of the total number of macroscopic grid points)


&maxwell
--------

- **al_em(3)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Size of simulation box in electromagnetic analysis. Unit of the length can be chosen by ``&units/unit_system``.

- **dl_em(3)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Spacing of real-space grids in electromagnetic analysis. Unit of length can be chosen by ``&units/unit_system``.

- **dt_em** (real(8), Default=0)
   | Available for ``theory='maxwell'``.
   Time step in electromagnetic analysis. Unit of time can be chosen by ``&units/unit_system``.

- **nt_em** (integer, Default=0)
   | Available for ``theory='maxwell'``.
   Number of total time steps for real-time propagation in electromagnetic analysis.

- **boundary_em(3,2)** (character, Default='default')
   | Available for ``theory='maxwell'`` and ``theory='multi_scale_maxwell_tddft'``.
   Boundary condition in electromagnetic analysis. The first index(1-3 rows) corresponds to x, y, and z axes. The second index(1-2 columns) corresponds to bottom and top of the axes. If ``&system/yn_periodic='n'``, ``'default'``, ``'abc'``, and ``'pec'`` can be chosen, where ``'default'`` automatically chooses ``'abc'``. If ``&system/yn_periodic='y'``, ``'default'``, ``'pml'``, and ``'periodic'`` can be chosen, where ``'default'`` automatically chooses ``'periodic'``. ``'abc'`` is absorbing boundary, ``'pec'`` is perfect electric conductor, and ``'periodic'`` is periodic boundary.

- **shape_file** (character, Default='none')
   | Available for ``theory='maxwell'``.
   Name of shape file in electromagnetic analysis. The shape files can be generated by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **media_num** (integer, Default=0)
   | Available for ``theory='maxwell'``.
   Number of media in electromagnetic analysis.

- **media_type(:)** (character, Default='vacuum')
   | Available for ``theory='maxwell'``.
   Type of media in electromagnetic analysis. ``'vacuum'``, ``'constant media'``, ``'pec'``, and ``'lorentz-drude'`` can be chosen. If ``'lorentz-drude'`` is chosen, linear response calculation can be done by ``&emfield/ae_shape1 or ae_shape2='impulse'``.

- **epsilon_em(:)** (real(8), Default=1d0)
   | Available for ``theory='maxwell'``.
   Relative permittivity of the media in electromagnetic analysis.

- **mu_em(:)** (real(8), Default=1d0)
   | Available for ``theory='maxwell'``.
   Relative permeability of the media in electromagnetic analysis.

- **sigma_em(:)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Conductivity of the media in electromagnetic analysis.

- **pole_num_ld(:)** (integer, Default=1)
   | Available for ``theory='maxwell'``.
   Number of poles of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis.

- **omega_p_ld(:)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Plasma frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis.
   
- **f_ld(:,:)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Oscillator strength of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``media_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``.

- **gamma_ld(:,:)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Collision frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``media_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``.

- **omega_ld(:,:)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Oscillator frequency of the media for the case of ``type_media='lorentz-drude'`` in electromagnetic analysis. The first index is media id whose maximum value is determined by ``media_num``. The second index is pole id whose maximum value is determined by ``pole_num_ld``.

- **wave_input** (character, Default='none')
   | Available for ``theory='maxwell'``.
   If ``'source'``, the incident pulse in electromagnetic analysis is generated by the incident current source.

- **ek_dir1(3)/ek_dir2(3)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Propagation direction of the first/second pulse.

- **source_loc1(3)/source_loc2(3)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Location of the incident current source of the first/second pulse. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/yn_periodic='n'`` while ranges from ``0`` to ``al_em`` for ``&system/yn_periodic='y'``.

- **obs_num_em** (integer, Default=0)
   | Available for ``theory='maxwell'``.
   Number of observation point in electromagnetic analysis. From the obtained results, figure and animation files can be generated by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **obs_samp_em** (integer, Default=1)
   | Available for ``theory='maxwell'``.
   Sampling time-step of the observation in electromagnetic analysis.

- **obs_loc_em(:,3)** (real(8), Default=0d0)
   | Available for ``theory='maxwell'``.
   Location of the observation point in electromagnetic analysis. Note that the coordinate system ranges from ``-al_em/2`` to ``al_em/2`` for ``&system/yn_periodic='n'`` while ranges from ``0`` to ``al_em`` for ``&system/yn_periodic='y'``.

- **yn_obs_plane_em(:)** (character, Default='n')
   | Available for ``theory='maxwell'``.
   Enable(``'y'``)/disable(``'n'``). Output of the electrmagnetic fields on the planes (xy, yz, and xz planes) for each observation point. This option must be ``'y'`` for generating animation files by using SALMON utilities (https://salmon-tddft.jp/utilities.html).

- **yn_wf_em** (character, Default='y')
   | Available for ``theory='maxwell'``.
   Enable(``'y'``)/disable(``'n'``). Applying a window function for linear response calculation when ``&calculation/theory=maxwell``.

&analysis
---------

- **projection_option** (character, Default='no')
   | Available for TDDFT based options of ``theory``.
   | Methods of projection to analyze the excited states (e.g. the number of excited electrons.)
   | Options
   |   ``'no'`` / no projection.
   |   ``'gs'`` / projection to eigenstates of ground-state Hamiltonian.
   |   ``'rt'`` / projection to eigenstates of instantaneous Hamiltonian.

- **out_projection_step** (integer, Default=100)
   | Available for ``projection_option`` with TDDFT based options of ``theory``.
   | Old infomation: 3d XXX
   Interval time step of projection analysis 

- **nenergy** (integer, Default=1000)
   | Available for ``theory='XXX', 'maxwell'``.
   Number of energy grid points for frequency-domain analysis.
   This parameter is required when `'impulse'` is choosen in `&emfield/ae_shape1|2`.

- **de** (real(8), Default=0.01d0 eV)
   | Available for ``theory='XXX', 'maxwell'``.
   Energy grid size for analysis.
   This parameter is required when `'impulse'` is choosen in `&emfield/ae_shape1|2`.

- **yn_out_psi** (character, Default='n')
   | Available for ``theory='XXX'``.
   | Option for output of wavefunctions
   | Options
   |   ``'y'`` / enable.
   |   ``'n'`` / disable.
   For periodic system (``yn_periodic='y'``), it works only for ground state calculation. The converged wave functions of all orbitals with all k-points are printed in gs_wfn_cube or gs_wfn_vtk directory. The format is speficied by ``format3d``.  XXX need check XXX

- **yn_out_dos** (character, Default='n')
   | Available for ``theory='dft'``.
   | Option for output of density of state
   | Options
   |   ``'y'`` / enable.
   |   ``'n'`` / disable.

- **yn_out_dos_set_fe_origin** (character, Default='n')
   | Available for ``theory='XXX'``.
   If ``'y'``, the electron energy is shifted to fix the Fermi energy as zero point.
   For ``&system/yn_periodic`` is ``'n'``, `` out_dos_fshift`` is not used 
   if ``&system/nstate`` is equal to ``&system/nelec``/2.

- **out_dos_start** (real(8), Default=-1d10 eV)
   | Available for ``theory='XXX'``.
   Lower bound (energy) of the density of state spectra.
   If this value is lower than a specific value near the lowest energy level, 
   this value is overwritten by that value. 

- **out_dos_end** (real(8), Default=1d10 eV)
   | Available for ``theory='XXX'``.
   Upper bound (energy) of the density of state spectra.
   If this value is higher than a specific value near the highest energy level, 
   this value is overwritten by that value. 

- **out_dos_nenergy** (integer, Default=601)
   | Available for ``theory='xxx'``.
   Number of  energy points sampled in the density of state spectra.
 
- **out_dos_width** (real(8), Default=0.1d0 eV)
   | Available for ``theory='XXX'``.
   Smearing width used in the density of state spectra..

- **out_dos_function** (character, Default='gaussian')
   | Available for ``theory='XXX'``.
   Choise of smearing method for the density of state spectra..
   ``gaussian`` and ``lorentzian`` function are available.

- **yn_out_pdos** (character, Default='n')
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   If ``'y'``, projected density of state is output.

- **yn_out_dns** (character, Default='n')
   | Available for ``theory='XXX'``.
   If ``'y'``, the spatial electron density distribution at the ground state is output.

- **yn_out_dns_rt/out_dns_rt_step** (Character/Integer, Default='n')
   | Available for ``theory='XXX'``.
   If ``'y'``,  the spatiotemporal electron density distribution during real-time time-propagation is output every ``outdns_rt_step`` time steps.

- **yn_out_dns_ac_je/out_dns_ac_je_step** (Character/Integer, Default='n'/50)
   | Available for ``theory='single_scale_maxwell_tddft'``.
   If ``'y'``,  the electron density, vector potential, electronic current, and ionic coordinates are printed out every ``outdns_dns_ac_je_step`` time steps.
   The data written in binary format are divided to files corresponding to the space-grid parallelization number. 
  
- **yn_out_dns_trans/out_dns_trans_energy** (Character/Real(8), Default='n'/1.55d0eV)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 3d
   If ``'y'``, transition in different density from the ground state at specified field frequency omega(given by ``out_dns_trans_energy``) is calculated by drho(r,omega)=FT(rho(r,t)-rho_gs(r))/T.

- **yn_out_elf** (character, Default='n')
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   If ``'y'``, electron localization function is output.

- **yn_out_elf_rt/out_elf_rt_step** (Character/Integer,Default='n'/50)
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   If ``'y'``, electron localization function 
   during real-time time-propagation is output
   every ``out_elf_rt_step`` time steps.

- **yn_out_estatic_rt/out_estatic_rt_step** (Character/Integer, Default='n'/50)
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   If ``'y'``, static electric field
   during real-time time-propagation is output
   every ``out_estatic_rt_step`` time steps.

- **yn_out_rvf_rt/out_rvf_rt_step** (Character/Integer, Default='n'/10)
   | Available for TDDFT based options and 'dft_md' option of ``theory``.
   If ``'y'``, coordinates[A], velocities[au], forces[au] on atoms
   during time-propagation are printed in ``SYSname``\_trj.xyz every ``out_rvf_rt_step`` time steps.
   If ``yn_md='y'``, the printing option is automatically turned on.
   
- **yn_out_tm** (character, Default='n')[Trial]
   | Available for ``yn_periodic='y'`` with ``theory='dft'``.
   If ``'y'``, transition moments between occupied and virtual orbitals are printed into ``SYSname``\_tm.data after the ground state calculation.

- **out_ms_step** (integer, Default=100)
   | Available for ``theory='multi_scale_maxwell_tddft'``.
   Some information is printed every ``out_ms_step`` time step in the Maxwell + TDDFT multi-scale calculation.

- **format_voxel_data** (character, Default='cube')
   | Available for ``theory='XXX'``.
   File format for three-dimensional volumetric data.
   ``'avs'``, ``'cube'``, and ``'vtk'`` can be chosen.

- **nsplit_voxel_data** (integer, Default=1)
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   Number of separated files for three dimensional data.
   Effective only when ``format3d`` is ``'avs'``.
   ``numfiles_out_3d`` must be less than or equal to number of processes.


&poisson
--------

- **layout_multipole** (character, Default=3)
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   | A variable to determine how to put multipoles in the Hartree potential calculation.
   | Options:
   |  ``1``/ A single pole is put at the center.
   |  ``2``/ Multipoles are put at the center of atoms.
   |  ``3``/ Multipoles are put at the center of mass of electrons in prepared cuboids.

- **num_multipole_xyz(3)** (integer, Default=0)
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   Number of multipoles when this is ``3``. When default is set, number of multipoles is calculated automatically.

- **threshold_cg** (real(8), Default=1d-15 a.u.)
   | Available for ``theory='XXX'``.
   | XXX


&ewald
------

- **newald** (integer, Default=4)
   | Available for ``yn_periodic='y'`` with DFT/TDDFT based options of ``theory``.
   Parameter for Ewald method. 
   Short-range part of Ewald sum is calculated within ``newald`` th
   nearlist neighbor cells.

- **aewald** (real(8), Default=0.5d0)
   | Available for ``yn_periodic='y'`` with DFT/TDDFT based options of ``theory``.
   Square of range separation parameter for Ewald method in atomic unit. 

- **cutoff_r** (real(8), Default=-1d0)
   | Available for ``yn_periodic='y'`` with DFT/TDDFT based options of ``theory``.
   | XXX

- **cutoff_r_buff** (real(8), Default=2d0 a.u.)
   | Available for ``yn_periodic='y'`` with DFT/TDDFT based options of ``theory``.
   | XXX

- **cutoff_g** (real(8), Default=-1d0)
   | Available for ``yn_periodic='y'`` with DFT/TDDFT based options of ``theory``.
   | XXX

&opt[Trial]
-------------

- **nopt** (integer, Default=100)
   | Available for ``yn_opt='y'`` with ``theory='dft'``.
   The maximum step number of geometry optimization.

- **convrg_opt_fmax** (real(8), Default=1d-3)
   | Available for ``yn_opt='y'`` with ``theory='dft'``.
   Convergence threshold of geometry optimization in maximum force.

- **max_step_len_adjust** (real(8), Default=-1d0)
  | Available for ``yn_opt='y'`` with ``theory='dft'``.
  | XXX

  
&md[Trial]
-----------
- **ensemble** (character, Default='NVE')[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Ensemble in MD option: "NVE" or "NVT".

- **thermostat** (character, Default='nose-hoover')[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Thermostat in "NVT" option: (currently only ``nose-hoover``).

- **step_velocity_scaling** (integer, Default=-1)[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Time step interval for velocity-scaling. Velocity-scaling is applied if this is set to positive.

- **step_update_ps** (Integer, Default=10)[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Time step interval for updating pseudopotential (Larger number makes calculation time reduce greatly, but gets inaccurate) in case of ``yn_md=y``.

- **temperature0_ion_k** (real(8), Default=298.15d0)[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Setting temperature [K] for NVT ensemble, velocity scaling and generating initial velocities.

- **yn_set_ini_velocity** (character, Default='n')[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Initial velocities are set.

  - ``y``: Generate initial velocity with Maxwell-Bortzman distribution.
  - ``r``: Read initial velocity from file specified by keyword of ``file_ini_velocity``. This is, for example, used for restarting MD from the previous run. The last atomic coordinates and velocities are printed in ``SYSname``\_trj.xyz. (atomic coordinate also should be copied from the previous output and put in the next input file for restart)
    
- **file_ini_velocity** (character, Default='none')[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   File name for initial velocities. This is read when ``set_ini_velocity`` is ``'r'``. The format is simply vx(iatom) vy(iatom) vz(iatom) in each line. The order of atoms must be the same as the given coordinates in the main input file. In case of using nose-hoover thermostat, a thermostat variable should be put at the last line (all atomic unit). 

- **thermostat_tau** (real(8), Default=41.34d0 a.u. or 1d0 fs)[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Parameter in Nose-Hoover method: controlling time constant for temperature.
   Default is ``41.34[au] or 1.0[fs]``.

..
#XXX removed?#
- **seed_ini_velocity** (integer, Default=123)[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Random seed (integer number) to generate initial velocity if ``set_ini_velocity`` is set to y.
   Default is ``123``.
..

- **yn_stop_system_momt** (character, Default='n')[Trial]
   | Available for ``yn_md='y'`` or ``theory='dft_md'``.
   Center of mass is fixed every time step.


&code
-----

- **yn_want_stencil_hand_vectorization** (character, Default='y')
   | This option requests hand-vectorized optimization code of stencil in the hamiltonian calculation.
   | SALMON checks the calculation can be used the hand-vectorized code.
   | If failing it, SALMON will uses the typical implementation.

- **yn_want_communication_overlapping** (character, Default='n')
   | Available for ``theory='tddft*' or '*maxwell_tddft'``
   | This option requests computation/communication overlap algorithm to improve the performance of stencil in the hamiltonian calculation.
   | SALMON checks the calculation can be used the overlap algorithm.
   | If failing it, SALMON will uses the non-overlap algorithm.

- **stencil_openmp_mode** (character, Default='auto')
   | This option selects a OpenMP parallelization mode of stencil in the hamiltonian calculation.
   | ``auto``    / SALMON decides the parallelization target automatically.
   | ``orbital`` / OpenMP parallelization is applied to orbital (and k-point) loop.
   | ``rgrid``   / OpenMP parallelization is applied to real-space grid loop.
   
- **current_openmp_mode** (character, Default='auto')
   | This option selects a OpenMP parallelization mode of the current calculation.
   | ``auto``    / SALMON decides the parallelization target automatically.
   | ``orbital`` / OpenMP parallelization is applied to orbital (and k-point) loop.
   | ``rgrid``   / OpenMP parallelization is applied to real-space grid loop.

- **force_openmp_mode** (character, Default='auto')
   | This option selects a OpenMP parallelization mode of the force calculation.
   | ``auto``    / SALMON decides the parallelization target automatically.
   | ``orbital`` / OpenMP parallelization is applied to orbital (and k-point) loop.
   | ``rgrid``   / OpenMP parallelization is applied to real-space grid loop.


**Following variables are moved from the isolated part. Some of them may be added to common input, be combined to it, and be removed.**


&group_fundamental[Trial]
-------------------------
   
- **ntmg** (integer, Default=1)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   Number of multigrid calculation for gs. At the moment, there is a malfunction in this variable, and recovery is needed.

- **iwrite_projection** (integer, Default=0)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   A variable for projection. 

- **itwproj** (integer, Default=-1)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   The projection is calculated every ``itwproj`` time steps. 

- **iwrite_projnum** (integer, Default=0)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   There is a malfunction in this variable.

- **itcalc_ene** (integer, Default=10)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   Total energy is calculated every ``itcalc_ene`` time steps. There may be a malfunction in this variable.


&group_hartree[Trial]
----------------------

- **lmax_lmp** (integer, Default=4)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   A maximum angular momentum for multipole expansion in the Hartree-cg calculation. 


&group_others[Trial]
---------------------

- **num_projection** (Interger, Default=1)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   Number of orbitals for projections.

- **iwrite_projection_ob(200)** (Interger, Default=1, 2, 3, ..., 200)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
   Orbital number to be written as projections.

- **iwrite_projection_k(200)** (Interger, Default=1)[Trial]
   | Available for ``theory='XXX'``.
   | Old infomation: 0d
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
.. _&group_hartree: #group_hartree
.. _&group_others: #group_others



