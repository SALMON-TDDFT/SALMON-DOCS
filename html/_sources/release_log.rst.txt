###########################
Release History
###########################

Release Notes
-------------


* Release note of SALMON ver. 2:

  https://github.com/SALMON-TDDFT/SALMON2/releases

* Release note of SALMON ver. 0 and 1:

  https://github.com/SALMON-TDDFT/SALMON/releases


Details of Minor Changes
------------------------

The following is the history of fixed bugs and changes in models/inputs/outputs after releasing v.2.0.0. (not complete list currently)

Fixed bugs
==========

(Fixed in v.2.2.2)

* Standard output of the DFT part with ``convergence/=rho_dne`` was incorrect.
* Some external links, such as those to pseudopotential databases and LAPACK libraries, were incorrect.
* ``--arch=intel-avx512`` was incompatible with the latest versions of Intel compilers.
* Some testsuites (192_bulk_Si_pseudo_pspnc & 724_Si_diamond_bloch_ms) had issues related to compiler dependencies.
* Some ``.pspnc`` pseudopotentials with certain r-grids caused out-of-bounds errors. Index checking has been updated.

(Fixed in v.2.2.1)

* The OpenACC mode by newer versions of Nvidia HPC SDK was not supported.
* GNU Compiler Collection (GCC) was not supported.
* The sign of the polarizability in isolated systems (yn_periodic=n) was wrong.
* The definition of the excitation energy for the frozen Hamiltonian calculation (yn_fix_func=y) was wrong.
* File output of the atomic force for the multiscale mode (e.g. theory=multi_scale_maxwell_tddft and yn_out_rvf_rt=y) was wrong.
* File output of ``sysname_rt.data`` for very large systems was wrong.
* An array allocation for the projection_option mode was wrong.
* An array allocation for spin-noncollinear isolated systems (yn_periodic=n and spin=noncollinear) was wrong.

(Fixed in v.2.2.0)

* For the incident pulse for isolated systems (``yn_periodic=’n’``), the circular polarization was not defined and the CEP for Acos2 envelope was wrong.
* For spin-unpolarized systems (``spin=‘unpolarized’``), the initial value of the occupation rate was wrong when the electron number (``nelec``) is odd.
* When ``unit_system=‘A_eV_fs’`` is specified, ``temperature`` was mistakenly defined in the atomic unit.
* For isolated systems (``yn_periodic=’n’``), the multipole expansion for the boundary conditions of the Hartree potential (Poisson equation) was wrong. 
* For periodic systems (``yn_periodic=’y’``), ``al_vec[123]`` for orthogonal cells yielded unintended error. 
* For spin-noncollinear systems (``spin=‘noncollinear’``), the band-decomposition of the spin magnetization in ***_mag.data files was wrong.

(Fixed in v.2.1.0)

* Non-local term of the transition moment printed by the option of "yn_out_tm=y" has been fixed (just printing issue).
* Parallelization for orbitals for calculation of transition moment by "yn_out_tm=y" has been suported
* Bug of the segmentation fault occurred by "yn_ffte=y" with parallelization for orbitals using isolated system has been fixed.
* Reading of CIF file format for symmetry option has been improved. 
* Combination of non-uniform user-defined k-points and symmetry option has been supported.


(Fixed in v.2.0.2)

* Small noise on the total energy in TDDFT calculation (that is seen with weak pulse around e.g. I=1d9 W/cm2) has been removed.
* The printed absolute values of electron density in cube format has been fixed.. 
* Printing of the external field in TDDFT calculation of the isolated system has a bug in v.2.0.0 and v.2.0.1. It has been fixed in v.2.0.2.
* The file reading option of the external electric field in TDDFT calculation ("file_input1" in &emfield) has been fixed.
* Invalid occupation number printed in SYSNAME_ovlp.data file for projection option with non-uniform k-points has been fixed.

(Fixed in v.2.0.0)

* The imaginary part of wavefunction was not printed in cube format until v.1.2.1

(Fixed in v.?.?.?)

* Abnormal calculation that sometimes happens if zero value is included in the atomic coordinate in the input with "A_eV_fs" has been fiexd. 


Changes of models/inputs/outputs
================================

(v.2.2.2)

* Filenames of some output files are changed. 

  * ``dos.data`` --> ``<base_directory>/<sysname>_dos.data``
  * ``dns.cube`` --> ``<base_directory>/<sysname>_dns.cube``
  * etc.

* Some input keywords and options are added.

  * ``dk_shift``: Shift of the k-vector.
  * ``yn_out_rt_energy_components``: Switch for printing out the energy components such as the kinetic energy term (for TDDFT).
  * ``magdir_atom``: Initial values for the spin polarization (for DFT).
  * ``method_mixing='simple_potential'``: Simple mixing method for the local potential (for DFT).


(v.2.2.1)

* New input keyword for the preconditioning of CG method for accelerating the DFT computation 

  * ``yn_preconditioning``

* For the spin-noncollinear mode (spin='noncollinear'), some input keywords and output files are changed and added.

  * ``sysname_rt_spin.data``: the output file for the spin magnetization and spin current density.
  * Change the input keyword: ``out_magnetization_step`` --> ``out_rt_spin_step``.
  * New input keywords: ``yn_out_mag_decomposed_rt``, ``yn_out_spin_current_decomposed``, ``yn_out_spin_current_micro``.

* New input keyword to read a ``.cube`` file of the initial electron density for accelerating the DFT computation

  * ``method_init_density='read_dns_cube'``

(v.2.2.0)

* New theory options for SBE and Maxwell-SBE are added.

  * ``theory = 'sbe'``
  * ``theory = 'maxwell_sbe'``
  
* Input keywords for the Poisson equation solver are added.

  * ``method_poisson``
  * ``yn_fftw``
  
* New TDDFT analysis options are added. 

  * ``yn_fix_func``
  * ``projection_option=‘td’``
  * ``threshold_projection``
  * ``yn_out_intraband_current``
  * ``yn_out_current_decomposed``, ``out_current_decomposed_step``
  * ``yn_out_micro_je``, ``out_micro_je_step``
  

(v.2.1.0)

* Input variables for Spin-orbit coupling are added

  * "yn_spinorbit"
  * "spin = noncollinear"
  * "out_magnetization_step"

* New options for calculation of dielectric function and conductivity based on transition moments at the end of the GS calculation is added. The related input variables are

  * "yn_out_gs_sgm_eps"
  * "out_gs_sgm_eps_mu_nu"
  * "out_gs_sgm_eps_width"


(v.2.0.2)

* The definition of the total energy of the periodic system printed in TDDFT calculation has changed: 
  The electric field energy is included until v.2.0.1. It has not been included from v.2.0.2. 
* The directory names generated by "method_wf_distributor=slice" option have changed from v.2.0.2.
