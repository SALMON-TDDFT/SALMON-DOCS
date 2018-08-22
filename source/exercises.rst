
Exercises
====================

Getting started
---------------

Welcome to SALMON Exercises!

In these exercises, we explain the use of SALMON from the very
beginning, taking a few samples that cover applications of SALMON in
several directions. We assume that you are in the computational
environment of UNIX/Linux OS. First you need to download and install
SALMON in your computational environment. If you have not yet done it,
do it following the instruction, `download <download and `Install
and Run <Install_and_Run.

As described in `Install and Run <Install_and_Run, you are required
to prepare at least an input file and pseudopotential files to run
SALMON. In the following, we present input files for several sample
calculations and provide a brief explanation of the namelist variables
that appear in the input files. You may modify the input files to
execute for your own calculations. Pseudopotential files of elements
that appear in the samples are also attached. We also present
explanations of main output files.

We present 6 exercises.

First 3 exercises (Exercise-1 ~ 3) are for an isolated molecule,
acetylene C2H2. If you are interested in learning electron dynamics
calculations in isolated systems, please look into these exercises. In
SALMON, we usually calculate the ground state solution first. This is
illustrated in
`Exercise-1 <#Exercise-1:_Ground_state_of_C2H2_molecule. After
finishing the ground state calculation, two exercises of electron
dynamics calculations are prepared.
`Exercise-2 <#Exercise-2:_Polarizability_and_photoabsorption_of_C2H2_molecule
illustrates the calculation of linear optical responses in real time,
obtaining polarizability and photoabsorption of the molecule.
`Exercise-3 <#Exercise-3:_Electron_dynamics_in_C2H2_molecule_under_a_pulsed_electric_field
illustrates the calculation of electron dynamics in the molecule under a
pulsed electric field.

Next 2 exercises (Exercise-4 ~ 5) are for a crystalline solid, silicon.
If you are interested in learning electron dynamics calculations in
extended periodic systems, please look into these exercises. Since
ground state calculations of small unit-cell systems are not
computationally expensive and a time evolution calculation is usually
much more time-consuming than the ground state calculation, we recommend
to run the ground and the time evolution calculations as a single job.
The following two exercises are organized in that way.
`Exercise-4 <#Exercise-4:_Dielectric_function_of_crystalline_silicon
illustrates the calculation of linear response properties of crystalline
silicon to obtain the dielectric function.
`Exercise-5 <#Exercise-5:_Electron_dynamics_in_crystalline_silicon_under_a_pulsed_electric_field
illustrates the calculation of electron dynamics in the crystalline
silicon induced by a pulsed electric field.

The final exercise (Exercise-6) is for an irradiation and a propagation
of a pulsed light in a bulk silicon, coupling Maxwell equations for the
electromagnetic fields of the pulsed light and the electron dynamics in
the unit cells. This calculation is quite time-consuming and is
recommended to execute using massively parallel supercomputers.
`Exercise-6 <#Exercise-6:_Pulsed-light_propagation_through_a_silicon_thin_film
illustrates the calculation of a pulsed, linearly polarized light
irradiating normally on a surface of a bulk silicon.

C2H2 (isolated molecules)
-------------------------

Exercise-1: Ground state of C2H2 molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this exercise, we learn the calculation of the ground state solution
of acetylene (C2H2) molecule, solving the static Kohn-Sham equation.
This exercise will be useful to learn how to set up calculations in
SALMON for any isolated systems such as molecules and nanoparticles. It
should be noted that at present it is not possible to carry out the
geometry optimization in SALMON. Therefore, atomic positions of the
molecule are specified in the input file and are fixed during the
calculations.

Input files
^^^^^^^^^^^

To run the code, following files are used:

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *C2H2_gs.inp*                     | input file that contains namelist |
|                                   | variables and their values        |
+-----------------------------------+-----------------------------------+
| *C_rps.dat*                       | pseodupotential file for carbon   |
|                                   | atom                              |
+-----------------------------------+-----------------------------------+
| *H_rps.dat*                       | pseudopotential file for hydrogen |
|                                   | atom                              |
+-----------------------------------+-----------------------------------+

You may download the above 3 files (zipped file) from:

```Download``\ ````\ ``zipped``\ ````\ ``input``\ ````\ ``and``\ ````\ ``pseudopotential``\ ````\ ``files`` <media:C2H2_gs_input.zip

In the input file *C2H2_gs.inp*, namelists variables are specified. Most
of them are mandatory to execute the ground state calculation. We
present their explanations below:

```Explanations``\ ````\ ``of``\ ````\ ``input``\ ````\ ``files``\ ````\ ``(ground``\ ````\ ``state``\ ````\ ``of``\ ````\ ``C2H2``\ ````\ ``molecule)`` <Explanations_of_input_files_(ground_state_of_C2H2_molecule)

This will help you to prepare an input file for other systems that you
want to calculate. A complete list of the namelist variables that can be
used in the input file can be found in the downloaded file
*SALMON/manual/input_variables.md*.

Output files
^^^^^^^^^^^^

After the calculation, following output files are created in the
directory that you run the code,

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *C2H2_info.data*                  | information on ground state       |
|                                   | solution                          |
+-----------------------------------+-----------------------------------+
| *dns.cube*                        | a cube file for electron density  |
+-----------------------------------+-----------------------------------+
| *elf.cube*                        | electron localization function    |
|                                   | (ELF)                             |
+-----------------------------------+-----------------------------------+
| *psi1.cube*, *psi2.cube*, ...     | electron orbitals                 |
+-----------------------------------+-----------------------------------+
| *dos.data*                        | density of states                 |
+-----------------------------------+-----------------------------------+
| *pdos1.data*, *pdos2.data*, ...   | projected density of states       |
+-----------------------------------+-----------------------------------+
| *C2H2_gs.bin*                     | binary output file to be used in  |
|                                   | the real-time calculation         |
+-----------------------------------+-----------------------------------+

You may download the above files (zipped file, except for the binary
file *C2H2_gs.bin*) from:

https://salmon-tddft.jp/wiki/media:C2H2_gs_output.zip

Main results of the calculation such as orbital energies are included in
*C2H2_info.data*. Explanations of the *C2H2_info.data* and other output
files are described in:

```Explanations``\ ````\ ``of``\ ````\ ``output``\ ````\ ``files``\ ````\ ``(ground``\ ````\ ``state``\ ````\ ``of``\ ````\ ``C2H2``\ ````\ ``molecule)`` <Explanations_of_output_files_(ground_state_of_C2H2_molecule)


We show several image that are created from the output files.

+-----------------------------------+-----------------------------------+
| image                             | files used to create the image    |
+-----------------------------------+-----------------------------------+
| `highest occupied molecular       | *psi1.cube*, *psi2.cube*, ...     |
| orbital                           |                                   |
| (HOMO) <:File:HOMO.png#file   |                                   |
+-----------------------------------+-----------------------------------+
| `electron                         | *dns.cube*                        |
| density <:File:Dns.png#file   |                                   |
+-----------------------------------+-----------------------------------+
| `electron localization            | *elf.cube*                        |
| function <:File:Elf.png#file  |                                   |
+-----------------------------------+-----------------------------------+

Exercise-2: Polarizability and photoabsorption of C2H2 molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this exercise, we learn the linear response calculation in the
acetylene (C2H2) molecule, solving the time-dependent Kohn-Sham
equation. The linear response calculation provides the polarizability
and the oscillator strength distribution of the molecule. This exercise
should be carried out after finishing the ground state calculation that
was explained in
`Exercise-1 <#Exercise-1:_Ground_state_of_C2H2_molecule. In the
calculation, an impulsive perturbation is applied to all electrons in
the C2H2 molecule along the molecular axis which we take *z* axis. Then
a time evolution calculation is carried out without any external fields.
During the calculation, the electric dipole moment is monitored. After
the time evolution calculation, a time-frequency Fourier transformation
is carried out for the electric dipole moment to obtain the
frequency-dependent polarizability. The imaginary part of the
frequency-dependent polarizability is proportional to the oscillator
strength distribution and the photoabsorption cross section.

.. _input-files-1:

Input files
^^^^^^^^^^^

To run the code, the input file *C2H2_rt_response.inp* that contains
namelist variables and their values for the linear response calculation
is required. The binary file *C2H2_gs.bin* that is created in the ground
state calculation and pseudopotential files are also required. The
pseudopotential files should be the same as those used in the ground
state calculation.

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *C2H2_rt_response.inp*            | input file that contains namelist |
|                                   | variables and their values        |
+-----------------------------------+-----------------------------------+
| *C_rps.dat*                       | pseodupotential file for carbon   |
+-----------------------------------+-----------------------------------+
| *H_rps.dat*                       | pseudopotential file for hydrogen |
+-----------------------------------+-----------------------------------+
| *C2H2_gs.bin*                     | binary file created in the ground |
|                                   | state calculation                 |
+-----------------------------------+-----------------------------------+

You may download the *C2H2_rt_response.inp* file (zipped file) from:

```Download``\ ````\ ``zipped``\ ````\ ``input``\ ````\ ``file`` <media:C2H2_rt_response_input.zip

In the input file *C2H2_rt_response.inp*, namelists variables are
specified. Most of them are mandatory to execute the linear response
calculation. We present their explanations below:

```Explanations``\ ````\ ``of``\ ````\ ``input``\ ````\ ``files``\ ````\ ``(polarizability``\ ````\ ``and``\ ````\ ``photoabsorption``\ ````\ ``of``\ ````\ ``C2H2``\ ````\ ``molecule)`` <Explanations_of_input_files_(polarizability_and_photoabsorption_of_C2H2_molecule)

This will help you to prepare the input file for other systems that you
want to calculate. A complete list of the namelist variables that can be
used in the input file can be found in the downloaded file
*SALMON/manual/input_variables.md*.

.. _output-files-1:

Output files
^^^^^^^^^^^^

After the calculation, following output files are created in the
directory that you run the code,

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *C2H2_lr.data*                    | polarizability and oscillator     |
|                                   | strength distribution as          |
|                                   | functions of energy               |
+-----------------------------------+-----------------------------------+
| *C2H2_p.data*                     | components of dipole moment as    |
|                                   | functions of time                 |
+-----------------------------------+-----------------------------------+

You may download the above files (zipped file) from:

https://salmon-tddft.jp/wiki/media:C2H2_rt_response_output.zip

Explanations of the output files are given in:

```Explanations``\ ````\ ``of``\ ````\ ``output``\ ````\ ``files``\ ````\ ``(polarizability``\ ````\ ``and``\ ````\ ``photoabsorption``\ ````\ ``of``\ ````\ ``C2H2``\ ````\ ``molecule)`` <Explanations_of_output_files_(polarizability_and_photoabsorption_of_C2H2_molecule)

Exercise-3: Electron dynamics in C2H2 molecule under a pulsed electric field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this exercise, we learn the calculation of the electron dynamics in
the acetylene (C2H2) molecule under a pulsed electric field, solving the
time-dependent Kohn-Sham equation. As outputs of the calculation, such
quantities as the total energy and the electric dipole moment of the
system as functions of time are calculated. This tutorial should be
carried out after finishing the ground state calculation that was
explained in
`Exercise-1 <#Exercise-1:_Ground_state_of_C2H2_molecule. In the
calculation, a pulsed electric field that has cos^2 envelope shape is
applied. The parameters that characterize the pulsed field such as
magnitude, frequency, polarization direction, and carrier envelope phase
are specified in the input file.

.. _input-files-2:

Input files
^^^^^^^^^^^

To run the code, following files are used. The *C2H2_gs.bin* file is
created in the ground state calculation. Pseudopotential files are
already used in the ground state calculation. Therefore,
*C2H2_rt_pulse.inp* that specifies namelist variables and their values
for the pulsed electric field calculation is the only file that the
users need to prepare.

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *C2H2_rt_pulse.inp*               | input file that contain namelist  |
|                                   | variables and their values.       |
+-----------------------------------+-----------------------------------+
| *C_rps.dat*                       | pseodupotential file for Carbon   |
+-----------------------------------+-----------------------------------+
| *H_rps.dat*                       | pseudopotential file for Hydrogen |
+-----------------------------------+-----------------------------------+
| *C2H2_gs.bin*                     | binary file created in the ground |
|                                   | state calculation                 |
+-----------------------------------+-----------------------------------+

You may download the *C2H2_rt_pulse.inp* file (zipped file) from:

```Download``\ ````\ ``zipped``\ ````\ ``input``\ ````\ ``file`` <media:C2H2_rt_pulse_input.zip

In the input file *C2H2_rt_pulse.inp*, namelists variables are
specified. Most of them are mandatory to execute the calculation of
electron dynamics induced by a pulsed electric field. We present
explanations of the namelist variables that appear in the input file in:

```Explanations``\ ````\ ``of``\ ````\ ``input``\ ````\ ``files``\ ````\ ``(C2H2``\ ````\ ``molecule``\ ````\ ``under``\ ````\ ``a``\ ````\ ``pulsed``\ ````\ ``electric``\ ````\ ``field)`` <Explanations_of_input_files_(C2H2_molecule_under_a_pulsed_electric_field)

This will help you to prepare the input file for other systems and other
pulsed electric fields that you want to calculate. A complete list of
the namelist variables that can be used in the input file can be found
in the downloaded file *SALMON/manual/input_variables.md*.

.. _output-files-2:

Output files
^^^^^^^^^^^^

After the calculation, following output files are created in the
directory that you run the code,

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *C2H2_p.data*                     | components of the electric dipole |
|                                   | moment as functions of time       |
+-----------------------------------+-----------------------------------+
| *C2H2_ps.data*                    | power spectrum that is obtained   |
|                                   | by a time-frequency Fourier       |
|                                   | transformation of the electric    |
|                                   | dipole moment                     |
+-----------------------------------+-----------------------------------+

You may download the above files (zipped file) from:

https://salmon-tddft.jp/wiki/media:C2H2_rt_pulse_output.zip

Explanations of the files are described in:

```Explanations``\ ````\ ``of``\ ````\ ``output``\ ````\ ``files``\ ````\ ``(C2H2``\ ````\ ``molecule``\ ````\ ``under``\ ````\ ``a``\ ````\ ``pulsed``\ ````\ ``electric``\ ````\ ``field)`` <Explanations_of_output_files_(C2H2_molecule_under_a_pulsed_electric_field)

Crystalline silicon (periodic solids)
-------------------------------------

Exercise-4: Dielectric function of crystalline silicon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this exercise, we learn the linear response calculation of the
crystalline silicon of a diamond structure. Calculation is done in a
cubic unit cell that contains eight silicon atoms. Since the ground
state calculation costs much less computational time than the time
evolution calculation, both calculations are successively executed.
After finishing the ground state calculation, an impulsive perturbation
is applied to all electrons in the unit cell along *z* direction. Since
the dielectric function is isotropic in the diamond structure,
calculated dielectric function should not depend on the direction of the
perturbation. During the time evolution, electric current averaged over
the unit cell volume is calculated. A time-frequency Fourier
transformation of the electric current gives us a frequency-dependent
conductivity. The dielectric function may be obtained from the
conductivity using a standard relation.

.. _input-files-3:

Input files
^^^^^^^^^^^

To run the code, following files are used:

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *Si_gs_rt_response.inp*           | input file that contain namelist  |
|                                   | variables and their values.       |
+-----------------------------------+-----------------------------------+
| *Si_rps.dat*                      | pseodupotential file of silicon   |
+-----------------------------------+-----------------------------------+

You may download the above 2 files (zipped file) from:

```Download``\ ````\ ``zipped``\ ````\ ``input``\ ````\ ``and``\ ````\ ``pseudopotential``\ ````\ ``files`` <media:_Si_gs_rt_response_input.zip

In the input file *Si_gs_rt_response.inp*, namelists variables are
specified. Most of them are mandatory to execute the calculation. We
present explanations of the namelist variables that appear in the input
file in:

```Explanations``\ ````\ ``of``\ ````\ ``input``\ ````\ ``files``\ ````\ ``(dielectric``\ ````\ ``function``\ ````\ ``of``\ ````\ ``crystalline``\ ````\ ``silicon)`` <Explanations_of_input_files_(dielectric_function_of_crystalline_silicon)

This will help you to prepare the input file for other systems that you
want to calculate. A complete list of the namelist variables that can be
used in the input file can be found in the downloaded file
*SALMON/manual/input_variables.md*.

.. _output-files-3:

Output files
^^^^^^^^^^^^

After the calculation, following output files are created in the
directory that you run the code,

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *Si_gs_info.data*                 | information of ground state       |
|                                   | calculation                       |
+-----------------------------------+-----------------------------------+
| *Si_eigen.data*                   | energy eigenvalues of orbitals    |
+-----------------------------------+-----------------------------------+
| *Si_k.data*                       | information on k-points           |
+-----------------------------------+-----------------------------------+
| *Si_rt.data*                      | electric field, vector potential, |
|                                   | and current as functions of time  |
+-----------------------------------+-----------------------------------+
| *Si_force.data*                   | force acting on atoms             |
+-----------------------------------+-----------------------------------+
| *Si_lr.data*                      | Fourier spectra of the dielectric |
|                                   | functions                         |
+-----------------------------------+-----------------------------------+
| *Si_gs_rt_response.out*           | standard output file              |
+-----------------------------------+-----------------------------------+

You may download the above files (zipped file) from:

https://salmon-tddft.jp/wiki/media:Si_gs_rt_response_output.zip

Explanations of the output files are described in:

```Explanation``\ ````\ ``of``\ ````\ ``output``\ ````\ ``fiels``\ ````\ ``(dielectric``\ ````\ ``function``\ ````\ ``of``\ ````\ ``crystalline``\ ````\ ``silicon)`` <Explanation_of_output_fiels_(dielectric_function_of_crystalline_silicon)

Exercise-5: Electron dynamics in crystalline silicon under a pulsed electric field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this exercise, we learn the calculation of electron dynamics in a
unit cell of crystalline silicon of a diamond structure. Calculation is
done in a cubic unit cell that contains eight silicon atoms. Since the
ground state calculation costs much less computational time than the
time evolution calculation, both calculations are successively executed.
After finishing the ground state calculation, a pulsed electric field
that has cos^2 envelope shape is applied. The parameters that
characterize the pulsed field such as magnitude, frequency,
polarization, and carrier envelope phase are specified in the input
file.

.. _input-files-4:

Input files
^^^^^^^^^^^

To run the code, following files are used:

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *Si_gs_rt_pulse.inp*              | input file that contain namelist  |
|                                   | variables and their values.       |
+-----------------------------------+-----------------------------------+
| *Si_rps.dat*                      | pseodupotential file for Carbon   |
+-----------------------------------+-----------------------------------+

You may download the above 2 files (zipped file) from:

```Download``\ ````\ ``zipped``\ ````\ ``input``\ ````\ ``and``\ ````\ ``pseudopotential``\ ````\ ``files`` <media:Si_gs_rt_pulse_input.zip

In the input file *Si_gs_rt_pulse.inp*, namelists variables are
specified. Most of them are mandatory to execute the calculation. We
present explanations of the namelist variables that appear in the input
file in:

```Explanation``\ ````\ ``of``\ ````\ ``input``\ ````\ ``files``\ ````\ ``(crystalline``\ ````\ ``silicon``\ ````\ ``under``\ ````\ ``a``\ ````\ ``pulsed``\ ````\ ``electric``\ ````\ ``field)`` <Explanation_of_input_files_(crystalline_silicon_under_a_pulsed_electric_field)

This will help you to prepare the input file for other systems that you
want to calculate. A complete list of the namelist variables that can be
used in the input file can be found in the downloaded file
*SALMON/manual/input_variables.md*.

.. _output-files-4:

Output files
^^^^^^^^^^^^

After the calculation, following output files are created in the
directory that you run the code,

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *Si_gs_info.data*                 | information of ground state       |
|                                   | calculation                       |
+-----------------------------------+-----------------------------------+
| *Si_eigen.data*                   | energy eigenvalues of orbitals    |
+-----------------------------------+-----------------------------------+
| *Si_k.data*                       | information on k-points           |
+-----------------------------------+-----------------------------------+
| *Si_rt.data*                      | electric field, vector potential, |
|                                   | and current as functions of time  |
+-----------------------------------+-----------------------------------+
| *Si_force.data*                   | force acting on atoms             |
+-----------------------------------+-----------------------------------+
| *Si_lr.data*                      | Fourier transformations of        |
|                                   | various quantities                |
+-----------------------------------+-----------------------------------+
| *Si_gs_rt_pulse.out*              | standard output file              |
+-----------------------------------+-----------------------------------+

You may download the above files (zipped file) from:

https://salmon-tddft.jp/wiki/media:Si_gs_rt_pulse_output.zip

Explanations of the output files are described in:

```Explanation``\ ````\ ``of``\ ````\ ``output``\ ````\ ``files``\ ````\ ``(crystalline``\ ````\ ``silicon``\ ````\ ``under``\ ````\ ``a``\ ````\ ``pulsed``\ ````\ ``electric``\ ````\ ``field)`` <Explanation_of_output_files_(crystalline_silicon_under_a_pulsed_electric_field)

Maxwell + TDDFT multiscale simulation
-------------------------------------

Exercise-6: Pulsed-light propagation through a silicon thin film
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this exercise, we learn the calculation of the propagation of a
pulsed light through a thin film of crystalline silicon. We consider a
silicon thin film of 53 nm thickness, and an irradiation of a few-cycle,
linearly polarized pulsed light normally on the thin film. First, to set
up initial orbitals, the ground state calculation is carried out. The
pulsed light locates in the vacuum region in front of the thin film. The
parameters that characterize the pulsed light such as magnitude and
frequency are specified in the input file. The calculation ends when the
reflected and transmitted pulses reach the vacuum region.

.. _input-files-5:

Input files
^^^^^^^^^^^

To run the code, following files are used:

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *Si_gs_rt_multiscale.inp*         | input file that contain namelist  |
|                                   | variables and their values.       |
+-----------------------------------+-----------------------------------+
| *Si_rps.dat*                      | pseodupotential file for silicon  |
+-----------------------------------+-----------------------------------+

You may download the above two files (zipped file) from:

```Download``\ ````\ ``zipped``\ ````\ ``input``\ ````\ ``and``\ ````\ ``pseudopotential``\ ````\ ``files`` <media:_Si_gs_rt_multiscale_input.zip

In the input file *Si_gs_rt_multiscale.inp*, namelists variables are
specified. Most of them are mandatory to execute the calculation. We
present explanations of the namelist variables that appear in the input
file in:

```Explanation``\ ````\ ``of``\ ````\ ``input``\ ````\ ``files``\ ````\ ``(pulsed-light``\ ````\ ``propagation``\ ````\ ``through``\ ````\ ``a``\ ````\ ``silicon``\ ````\ ``thin``\ ````\ ``film)`` <Explanation_of_input_files_(pulsed-light_propagation_through_a_silicon_thin_film)

This will help you to prepare the input file for other systems that you
want to calculate. A complete list of the namelist variables that can be
used in the input file can be found in the downloaded file
*SALMON/manual/input_variables.md*.

.. _output-files-5:

Output files
^^^^^^^^^^^^

After the calculation, new directory *multiscale/* is created, then,
following output files are created in the directory,

+-----------------------------------+-----------------------------------+
| file name                         | description                       |
+-----------------------------------+-----------------------------------+
| *Si_gs_info.data*                 | results of the ground state as    |
|                                   | well as input parameters          |
+-----------------------------------+-----------------------------------+
| *Si_eigen.data*                   | orbital energies in the ground    |
|                                   | state calculation                 |
+-----------------------------------+-----------------------------------+
| *Si_k.data*                       | information on k-points           |
+-----------------------------------+-----------------------------------+
| *RT_Ac/Si_Ac_xxxxxx.data*         | various quantities at a time as   |
|                                   | functions of macroscopic position |
+-----------------------------------+-----------------------------------+
| *RT_Ac/Si_Ac_vac.data*            | vector potential at vacuum        |
|                                   | position adjacent to the medium   |
+-----------------------------------+-----------------------------------+
| *Mxxxxxx/Si_Ac_M.data*            | various quantities at a           |
|                                   | macroscopic point as functions of |
|                                   | time                              |
+-----------------------------------+-----------------------------------+
| *Si_gs_rt_multiscale.out*         | standard output file              |
+-----------------------------------+-----------------------------------+

You may download the above files (zipped file) from:

https://salmon-tddft.jp/wiki/media:Si_gs_rt_multiscale_output.zip

Explanations of the output files are described in:

```Explanation``\ ````\ ``of``\ ````\ ``output``\ ````\ ``files``\ ````\ ``(pulsed-light``\ ````\ ``propagation``\ ````\ ``through``\ ````\ ``a``\ ````\ ``silicon``\ ````\ ``thin``\ ````\ ``film)`` <Explanation_of_output_files_(pulsed-light_propagation_through_a_silicon_thin_film)
