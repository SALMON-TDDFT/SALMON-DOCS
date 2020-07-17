###########################
Introduction
###########################

Overview
-----------

SALMON is an open-source computer program for ab-initio
quantum-mechanical calculations of electron dynamics at the nanoscale
that takes place in various situations of light-matter interactions. It
is based on time-dependent density functional theory, solving
time-dependent Kohn-Sham equation in real time and real space with
norm-conserving pseudopotentials.

SALMON was born by unifying two scientific programs: ARTED, developed by
Univ. Tsukuba group, that describes electron dynamics in crystalline
solids, and GCEED, developed by Institute for Molecular Science group,
that describes electron dynamics in molecules and nanostructures. It can
thus describe electron dynamics in both isolated and periodic systems.
It can also describe coupled dynamics of electrons and light-wave
electromagnetic fields.

To run the program, SALMON requires MPI Fortran/C compiller with LAPACK
libraries. SALMON has been tested and optimized to run in a number of
platforms, including Linux PC Cluster with x86-64 CPU, supercomputer
systems with Fujitsu FX100 and A64FX processors, and supercomputer system 
with Intel Xeon Phi (Knights Landing).

SALMON features
-------------------

In the microscopic scale, SALMON describes electron dynamics in both 
isolated (molecules and nanostructures) and periodic (crystalline solids) 
systems, solving time-dependent Kohn-Sham equation in real time and real space
with norm-conserving pseudopotential.
SALMON first carries out ground-state calculations in the density functional theory
to prepare initial configurations. SALMON then calculates electron
dynamics induced by applied electric field. Employing a weak impulsive
external field, SALMON can be used to calculate linear response
properties such as a polarizability of molecules and a dielectric
function of crystalline solids. Using pulsed electric fields, SALMON
describes electron dynamics in matters induced by intense and ultrashort
laser pulses.

SALMON is also capable of describing a propagation of electromagnetic fields 
of light using finite-difference time-domain method. As a unique feature
of SALMON, it is possible to carry out calculations of a coupled dynamics
of light electromagnetic fields and electron dynamics simultaneously.


Efficient parallelizations are implemented in the code by dividing spatial
grids, orbital index, and k-points. 
SALMON shows a good scalability when it runs in parallel supercomputers,
both for the ground state and the time evolution calculations.

-  Ground state calculations

   -  Kohn-Sham orbitals and energies
   -  density of states
   -  projected density of states
   -  electron localization function

-  Optical properties

   -  Oscillator strength distribution (absorption spectrum)
   -  dielectric function

-  Light-induced electron dynamics

   -  time evolution of Kohn-Sham orbitals
   -  density, current
   -  excitation energy
   -  number density of excited carriers

-  Propagation of light electromagnetic fields

   - Drude-Lorentz model
   - optical response of metasurfaces

-  Simultaneous description of electron dynamics and light pulse
   propagation

   -  light pulse propagation as well as time evolution of Kohn-Sham
      orbitals
   -  energy transfer from pulsed light to electrons


License
-----------

SALMON is available under Apache License version 2.0.

  Copyright 2017 SALMON developers

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at 

  http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and 
  limitations under the License.


  
SALMON at Github
--------------------

SALMON is developed at `GitHub.com <https://github.com/salmon-tddft>`__

List of developers
----------------------

(Alphabetic order)


* Yuta Hirokawa (University of Tsukuba, Japan)
* Kenji Iida (Hokkaido University, Japan)
* Tomohito Otobe (National Institutes for Quantum and Radiological Science and Technology, Japan)
* Shunsuke Sato (University of Tsukuba, Japan)
* Yasushi Shinohara (University of Tokyo, Japan)
* Takashi Takeuchi (University of Tsukuba, Japan)
* Mitsuharu Uemoto (Kobe University, Japan)
* Kazuhiro Yabana (University of Tsukuba, Japan)
* Atsushi Yamada (University of Tsukuba, Japan)
* Shunsuke Yamada (University of Tsukuba, Japan)

Former developers
----------------------

* Isabella Floss
* Kazuya Ishimura
* Kyung-Min Lee
* Katsuyuki Nobusada
* Masashi Noda
* Xiao-Min Tong
* Maiku Yamaguchi

Acknowledgements for SALMON developments
--------------------------------------------

SALMON has been developed by the SALMON developers under supports by
Center for Computational Sciences, University of Tsukuba, and 
National Institute for Quantum and Radiological Science and Technology.
SALMON has been supported by Strategic Basic
Research Programs, CREST, Japan Science and Technology Agency, under the
Grand Number JPMJCR16N5, in the research area of Advanced core
technology for creation and practical utilization of innovative
properties and functions based upon optics and photonics. SALMON was
also supported by Ministry of Education, Culture, Sports and
Technology of Japan as a social and scientific priority issue (Creation
of new functional devices and high-performance materials to support
next-generation industries: CDMSI) to be tackled by using post-K
computer.

