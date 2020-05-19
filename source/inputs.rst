.. _Inputs:

Inputs
====================

We here summarize input keywords that appear in :any:`Exercise <Exercises>`. 
A thorough list of the input keywords can be found in the downloaded file in
:any:`List of all input keywords <List of all input keywords>`.

.. _&calculation:

&calculation
------------

Mandatory: theory

::

   &calculation
     theory = 'dft'
   /

The value of the ``theory`` should be one of ``'dft'``, ``'dft_md'``, ``'tddft_response'``, ``'tddft_pulse'``,
``'single_scale_maxwell_tddft'``, ``'multi_scale_maxwell_tddft'``, ``'maxwell'``, and ``'dft_k_expand'``.

``dft`` : ground state calculation based on DFT

``dft_md`` : adiabatic ab initio MD simulations based on DFT

``tddft_response`` : simulations based on TDDFT to obtain linear-response

``tddft_pulse`` : simulations based on TDDFT using pulsed light

``single_scale_maxwell_tddft`` : coupled Maxwell and TDDFT single-scale simulation

``multi_scale_maxwell_tddft`` : coupled Maxwell and TDDFT multi-scale simulation

``maxwell`` : electromagnetic simulations based on the Maxwell's equations

``dft_k_expand`` : convert checkpoint data of dft with k-points calculation to that of larger supercell system with gamma-point

.. _&units:

&units
------

Mandatory: none

::

   &units
     unit_system = 'A_eV_fs'
   /

This namelist specifies the unit system to be used in the input file.
Options are 'A_eV_fs' for Angstrom, eV, and fs, and 'a.u.' or 'au' for atomic units.
If you do not specify it, atomic unit will be used as default.

For isolated systems (specified by ``yn_periodic = 'n'`` in ``&system``),
the unit of 1/eV is used for the output files of DOS and PDOS
if ``unit_system = 'A_eV_fs'`` is specified, while atomic unit is used if not. 
For other output files, the Angstrom/eV/fs units are used irrespective of the input keyword.

For periodic systems (specified by ``yn_periodic = 'y'`` in ``&system``),
the unit system specified by this input keyword is used for most output files.
See the first few lines of output files to confirm the unit system adopted in the file.

.. _&control:

&control
--------

Mandatory: none

::

   &control
     sysname = 'C2H2'
   /

'C2H2' defined by ``sysname = 'C2H2'`` will be used in the filenames of output files.
If you do not specify it, the file name will start with 'default'.

.. _&system:

&system
-------

Mandatory: yn_periodic, al, nelem, natom, nelec, nstate 

**For an isolated molecule (Exercises-1, 2, 3, 8, 9)**:

::

   &system
     yn_periodic = 'n'
     al(1:3) = 16.0d0, 16.0d0, 16.0d0
     nelem = 2
     natom = 4
     nelec = 10
     nstate = 6
   /

``yn_periodic = 'n'`` indicates that the isolated boundary condition will be used in the calculation.
``al(1:3) = 16.0d0, 16.0d0, 16.0d0`` specifies the lengths of three sides of the rectangular parallelepiped
where the grid points are prepared.
``nelem = 2`` and ``natom = 4`` indicate the number of elements and the number of atoms in the system, respectively.
``nelec = 10`` indicate the number of valence electrons in the system.
 ``nstate = 6`` indicates the number of Kohn-Sham orbitals to be solved.
``nstate`` should be equal to or larger than ``nelec/2``.

**For a periodic system (Exercises-4, 5, 6, 7)**:

::

   &system
     yn_periodic = 'y'
     al(1:3) = 10.26d0, 10.26d0, 10.26d0
     nelem = 1
     natom = 8
     nelec = 32
     nstate = 32
   /

``yn_periodic = 'y'`` indicates that three dimensional periodic boundary condition (bulk crystal) is assumed.
``al(1:3) = 10.26d0, 10.26d0, 10.26d0`` specifies the lattice constans of the unit cell.
``nelem = 1`` and ``natom = 8`` indicate the number of elements and the number of atoms in the system, respectively.
``nelec = 32`` indicate the number of valence electrons in the system. 
``nstate = 32`` indicates the number of Kohn-Sham orbitals to be solved.

.. _&pseudo:

&pseudo
-------

Mandatory: pseudo_file, izatom

**For C2H2 molecule**:

::

   &pseudo
     file_pseudo(1) = './C_rps.dat'
     file_pseudo(2) = './H_rps.dat'
     izatom(1) = 6
     izatom(2) = 1
     lloc_ps(1) = 1
     lloc_ps(2) = 0
   /

Parameters related to atomic species and pseudopotentials.
``pseudo_file(1) = './C_rps.dat'`` indicates the filename of the pseudopotential of element.
``izatom(1) = 6`` specifies the atomic number of the element.

**For crystalline Si**:

::

   &pseudo
     file_pseudo(1) = './Si_rps.dat'
     izatom(1) = 14
     lloc_ps(1) = 2
   /

``file_pseudo(1) = './Si_rps.dat'`` indicates the pseudopotential filename of element.
``izatom(1) = 14`` indicates the atomic number of the element.

.. _&functional:

&functional
-----------

::

   &functional
     xc ='PZ'
   /

``xc ='PZ'`` indicates that (adiabatic) local density approximation is adopted (Perdew-Zunger: Phys. Rev. B23, 5048 (1981)).
This is the default choice.

For isolated systems (specified by ``yn_periodic = 'n'`` in ``&system``),
only the default choice of 'PZ' is available at present.

For periodic systems (specified by ````yn_periodic = 'y'`` in ``&system``),
the following functionals may be available in addition to 'PZ', ``xc = 'PZM'``

Perdew-Zunger LDA with modification to improve sooth connection between
high density form and low density one, ``xc = 'TBmBJ' cval = 1.0``
:J. P. Perdew and Alex Zunger, Phys. Rev. B 23, 5048 (1981).

Tran-Blaha meta-GGA exchange with Perdew-Wang correlation. :Fabien Tran
and Peter Blaha, Phys. Rev. Lett. 102, 226401 (2009). John P. Perdew and
Yue Wang, Phys. Rev. B 45, 13244 (1992). This potential is known to
provide a reasonable description for the bandage of various insulators.
For this choice, the additional mixing parameter 'cval' may be
specified. If cval is set to a minus value, the mixing-parameter will be
computed following the formula in the original paper [Phys. Rev. Lett.
102, 226401 (2009)]. The default value for this parameter is 1.0.

Since version 1.1.0, exchange-correlation functionals in Libxc library
(http://www.tddft.org/programs/libxc/) have been usable in SALMON. At
present, usable functionals are limited to LDA and GGA. For periodic
systems, meta-GGA functionals are usable as well. To specify the
exchange-correlation potentials of Libxc library, there are two ways. If
the exchange and correlation potentials are given separately, you need
to specify both ``alibx`` and ``alibc`` separately. If the exchange and
correlation potentials are given as a combined set, you need to specify
``alibxc``. We show below an example:

::

   &functional
     alibx = 'LDA_X'
     alibc = 'LDA_C_PZ'
   /

Available sets of the functionals are listed at the website
http://www.tddft.org/programs/libxc/functionals/ .

Note that, the hybrid functionals (hybrid gga/mgga) are not supported in the current (version 2.0.0) of SALMON.

.. _&rgrid:

&rgrid
------

Mandatory: dl or num_rgrid

``dl(1:3) = 0.25d0, 0.25d0, 0.25d0`` specify the grid spacing in three Cartesian coordinates.
This is adopted for C2H2 calculation (Excercises-1, 2, 3, 8, 9).

::

   &rgrid
     dl(1:3) = 0.25d0, 0.25d0, 0.25d0
   /

``num_rgrid(1:3) = 12, 12, 12`` specify the number of grid points in each Cartesian direction.
This is adopted for crystalline Is calculation (Excercises-4, 5, 6, 7).

::

   &rgrid
     num_rgrid(1:3) = 12, 12, 12
   /

.. _&kgrid:
   
&kgrid
------

Mandatory: none

This group provides grid spacing of k-space for periodic systems.

::

   &kgrid
     num_kgrid(1:3) = 4, 4, 4
   /

.. _&scf:
   
&scf
----

Mandatory: nscf

This group specifies parameters related to the self-consistent field calculation.

::

   &scf
     nscf = 200
     threshold = 1.0d-9
   /

``nscf = 200`` is the number of scf iterations in the ground state calculation.
the scf loop in the ground state calculation ends before the number of the scf iterations reaches ``nscf``,
if a convergence criterion is satisfied. 

.. _&analysis:

&analysis
---------

Mandatory: none

The following input keywords specify whether the output files are created or not after the calculation.
In the ground state calculation of isolated systems (Excercise-1):

::

   &analysis
     yn_out_psi  = 'y'
     yn_out_dns  = 'y'
     yn_out_dos  = 'y'
     yn_out_pdos = 'y'
     yn_out_elf  = 'y'
   /

In the following input keywords, variables related to time-frequency Fourier analysis are specified.

::

   &analysis
     de = 1.0d-2
     nenergy = 3000
   /

``de = 1.0d-2`` specifies the energy spacing, 
and ``nenergy = 3000`` specifies the number of energy steps
in the time-frequency Fourier transformation.

.. _&tgrid:

&tgrid
------

Mandatory: dt, nt

::

   &tgrid
     dt = 1.25d-3
     nt = 5000
   /

``dt = 1.25d-3`` specifies the time step of the time evolution calculation. 
``nt = 5000`` specifies the number of time steps in the calculation.

.. _&emfield:

&emfield
--------

This group specifies the pulse shape of an electric filed applied to
the system in time evolution calculations. We explain below separating
two cases, :any:`linear-response-calculations`
and :any:`pulsed-electric-field-calculations`.

.. _linear-response-calculations:

Linear response calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A weak impulsive field is applied at *t=0*. For this case,
``ae_shape1 = 'impulse'`` should be described.

Mandatory: ae_shape1

::

   &emfield
     ae_shape1 = 'impulse'
     epdir_re1 = 0.d0,0.d0,1.d0
   /

``epdir_rel(3)`` specify a unit vector that indicates the direction of
the impulse.

For a periodic system specified by ``iperiodic = 3``, one may add
``trans_longi``. It has the value, ``'tr'``\ (transverse) or
``'lo'``\ (longitudinal), that specifies the treatment of the
polarization in the time evolution calculation. The default is ``'tr'``.

::

   &emfield
     trans_longi = 'tr'
     ae_shape1 = 'impulse'
     epdir_re1 = 0.,0.,1.
   /

The magnitude of the impulse of the pulse may be explicitly specified
by, for example, ``e_impulse = 1d-2``. The default is '1d-2' in atomic
unit.

.. _pulsed-electric-field-calculations:

Pulsed electric field calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Pulsed electric field of finite time duration is applied. For this
case, ``as_shape1`` should be specified. It indicates the shape of the
envelope of the pulse. The options include 'Acos2' and 'Ecos2' (See
below for other options).

Mandatory: ae_shape1, epdir_re1, {rlaser_int1 or amplitude1}, omega1,
pulse_tw1, phi_cep1

::

   &emfield
     ae_shape1 = 'Ecos2'
     epdir_re1 = 0.d0,0.d0,1.d0
     rlaser_int_wcm2_1 = 1.d8
     omega1=9.28d0
     pulse_tw1=6.d0
     phi_cep1=0.75d0
   /

``ae_shape1 = 'Ecos2'`` specifies the envelope of the pulsed electric
field, 'Ecos2' for the cos^2 envelope for the electric field. If 'Acos2'
is specified, this gives cos^2 envelope for the vector potential. Note
that 'phi_cep1' must be 0.75 (or 0.25) if one employs 'Ecos2' pulse
shape, since otherwise the time integral of the electric field does not
vanish. There is no such restriction for the 'Acos2' pulse shape.

``epdir_re1 = 0.d0,0.d0,1.d0`` specifies the real part of the unit
polarization vector of the pulsed electric field. If only the real part
is specified, it describes a linearly polarized pulse. Using both real
('epdir_re1') and imaginary ('epdir_im1') parts of the polarization
vector, circularly (and general ellipsoidary) polarized pulses may be
described.

``laser_int_wcm2_1 = 1.d8`` specifies the maximum intensity of the
applied electric field in unit of W/cm^2. It is also possible to specify
the maximum intensity of the pulse by ``amplitude1``.

``omega1=9.26d0`` specifies the average photon energy (frequency
multiplied with hbar).

``pulse_tw1=6.d0`` specifies the pulse duration. Note that it is not the
FWHM but a full duration of the cos^2 envelope.

``phi_cep1=0.75d0`` specifies the carrier envelope phase of the pulse.
As noted above, 'phi_cep1' must be 0.75 (or 0.25) if one employs 'Ecos2'
pulse shape, since otherwise the time integral of the electric field
does not vanish. There is no such restriction for the 'Acos2' pulse
shape.

It is possible to use two pulses simultaneously to simulate pump-probe
experiments, adding information for two pulses. To specify the second
pulse, change from 1 to 2 in the namelist variables, like ``ae_shape2``.
The time delay between two pulses is specified by the variable 't1_t2'.

For a periodic system specified by ``iperiodic = 3``, one may add
``trans_longi``. It has the value, ``'tr'``\ (transverse) or
``'lo'``\ (longitudinal), that specifies the treatment of the
polarization in the time evolution calculation. The default is ``'tr'``.
For a periodic system, it is also specify 'Acos3', 'Acos4', 'Acos6',
'Acos8' for ``ae_shape1``.

.. _&propagation:

&propagation
------------

This namelist specifies the numerical method for time evolution
calculations of electron orbitals.

::

   &propagation
     propagator='etrs'
   /

``propagator = 'etrs'`` indicates the use of enforced time-reversal
symmetry propagator. `M.A.L. Marques, A. Castro, G.F. Bertsch, and A.
Rubio, Comput. Phys. Commun., 151 60
(2003) <https://doi.org/10.1016/S0010-4655(02)00686-0>`__.

::

   &propagation
     propagator='middlepoint'
   /

``propagation='middlepoint'`` indicates that Hamiltonian at midpoint of
two-times is used.

The default is *middlepoint*.

.. _&multiscale:

&multiscale
-----------

This namelist specifies information necessary for Maxwell - TDDFT
multiscale calculations.

::

   &multiscale
     fdtddim = '1D'
     twod_shape = 'periodic'
     nx_m = 4
     ny_m = 1
     hX_m = 250d0
     nxvacl_m = -2000
     nxvacr_m = 256
   /

``fdtddim`` specifies the spatial dimension of the macro system.
``fdtddim='1D'`` indicates that one-dimensional equation is solved for
the macroscopic vector potential.

``nx_m = 4`` specifies the number of the macroscopic grid points in for
x-direction in the spatial region where the material exists.

``hx_m = 250d0`` specifies the grid spacing of the macroscopic grid in
x-direction.

``nxvacl_m = -2000`` and ``nxvacr_m = 256`` indicate the number of grid
points in the vacuum region, ``nxvacl_m`` for the left and ``nxvacr_m``
for the right from the surface of the material.

.. _&atomic_coor:

&atomic_coor
------------

Mandatory: atomic_coor or atomic_red_coor (they may be provided as a
separate file)

**For C2H2 molecule**:

::

   &atomic_coor
     'C' 0.000000 0.000000 0.599672 1
     'H' 0.000000 0.000000 1.662257 2
     'C' 0.000000 0.000000 -0.599672 1
     'H' 0.000000 0.000000 -1.662257 2
   /

Cartesian coordinates of atoms. The first column indicates the element.
Next three columns specify Cartesian coordinates of the atoms. The
number in the last column labels the element.


.. _&atomic_red_coor:

&atomic_red_coor
----------------

Mandatory: atomic_coor or atomic_red_coor (they may be provided as a
separate file)

**For a crystalline silicon**:

::

   &atomic_red_coor
     'Si' .0 .0 .0 1
     'Si' .25 .25 .25 1
     'Si' .5 .0 .5 1
     'Si' .0 .5 .5 1
     'Si' .5 .5 .0 1
     'Si' .75 .25 .75 1
     'Si' .25 .75 .75 1
     'Si' .75 .75 .25 1
   /

Cartesian coordinates of atoms are specified in a reduced coordinate
system. First column indicates the element, next three columns specify
reduced Cartesian coordinates of the atoms, and the last column labels
the element.



