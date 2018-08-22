###########################
Introduction
###########################

About SALMON
----------------

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
platforms, including Linux PC Cluster with x86-64 CPU, Fujitsu FX100
supercomputer system, K-computer, and supercomputer system with Intel
Xeon Phi (Knights Landing).

SALMON features
-------------------

SALMON describes electron dynamics in both isolated (molecules and
nanostructures) and periodic (crystalline solids) systems. SALMON first
carries out ground-state calculations in the density functional theory
to prepare initial configurations. SALMON then calculates electron
dynamics induced by applied electric field. Employing a weak impulsive
external field, SALMON can be used to calculate linear response
properties such as a polarizability of molecules and a dielectric
function of crystalline solids. Using pulsed electric fields, SALMON
describes electron dynamics in matters induced by intense and ultrashort
laser pulses.

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

-  Simultaneous description of electron dynamics and light pulse
   propagation

   -  light pulse propagation as well as time evolution of Kohn-Sham
      orbitals
   -  energy transfer from pulsed light to electrons

License
-----------

SALMON is available under Apache License version 2.0.

Copyright 2017 SALMON developers

Licensed under the Apache License, Version 2.0 (the "License"); you may
not use this file except in compliance with the License. You may obtain
a copy of the License at

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

* Isabella Floss (TU Wien, Austria)
* Yuta Hirokawa (University of Tsukuba, Japan)
* Kenji Iida (Institute for Molecular Science, Japan)
* Kazuya Ishimura (Institute for Molecular Science, Japan)
* Kyung-Min Lee (Max Planck Institute for the Structure and Dynamics of Matter, Germany)
* Katsuyuki Nobusada (Institute for Molecular Science, Japan)
* Masashi Noda (University of Tsukuba, Japan)
* Tomohito Otobe (National Institutes for Quantum and Radiological Science and Technology, Japan)
* Shunsuke Sato (Max Planck Institute for the Structure and Dynamics of Matter, Germany)
* Yasushi Shinohara (University of Tokyo, Japan)
* Takashi Takeuchi (University of Tsukuba, Japan)
* Xiao-Min Tong (University of Tsukuba, Japan)
* Mitsuharu Uemoto (University of Tsukuba, Japan)
* Kazuhiro Yabana (University of Tsukuba, Japan)
* Atsushi Yamada (University of Tsukuba, Japan)
* Shunsuke Yamada (University of Tsukuba, Japan)
* Maiku Yamaguchi (University of Tokyo, Japan)

Acknowledgements for SALMON developments
--------------------------------------------

SALMON has been developed by the SALMON developers under supports by
Center for Computational Sciences, University of Tsukuba, and Institute
for Molecular Science. SALMON has been supported by Strategic Basic
Research Programs, CREST, Japan Science and Technology Agency, under the
Grand Number JPMJCR16N5, in the research area of Advanced core
technology for creation and practical utilization of innovative
properties and functions based upon optics and photonics. SALMON has
been also supported by Ministry of Education, Culture, Sports and
Technology of Japan as a social and scientific priority issue (Creation
of new functional devices and high-performance materials to support
next-generation industries: CDMSI) to be tackled by using post-K
computer.
