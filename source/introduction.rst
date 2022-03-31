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

* Isabella Floss (TU Wien, Austria)
* Yuta Hirokawa (Prometech Software, Inc., Japan)
* Kenji Iida (Hokkaido University, Japan)
* Jun-Ichi Iwata (Advance Soft Co., Japan)
* Yuki Ito (Fixstars Corporation, Japan)
* Masashi Noda (Academeia, Japan)
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

* Kazuya Ishimura
* Kyung-Min Lee
* Katsuyuki Nobusada
* Xiao-Min Tong
* Maiku Yamaguchi


..
  We use sphinxcontrib-bibtex package for citing papers
  https://sphinxcontrib-bibtex.readthedocs.io/en/latest/index.html


.. _reference:

How to cite SALMON
--------------------

Suggested Citations
~~~~~~~~~~~~~~~~~~~~~
If you publish a paper in which SALMON makes an important contribution, please cite the SALMON code paper,
Ref. :cite:`Noda2019`
published in Computer Physics Communications.

We also suggest you to cite the following papers depending on your usage of SALMON.

* If you use SALMON for electron dynamics calculations of a large-size system,
  Ref. :cite:`Noda2014`
  that discusses massively parallel implementation utilizing spatial divisions will be appropriate.

* if you use SALMON to calculate electron dynamics in a unit cell of crystalline solid,
  Ref. :cite:`Bertsch2000`
  discussing formalism and numerical implementation will be appropriate.

* Ref. :cite:`Yabana1996`
  is one of the first implementations of the real-time time-dependent density functional calculation, in particular, instantaneous kick for the linear response calculations.

* If you use multiscale calculation coupling Maxwell equations for the electromagnetic fields of light and electron dynamics,
  Ref. :cite:`Yabana2012`
  discussing the formalism and the numerical implementation will be appropriate.

* Ref. :cite:`Sato2014JASSE`
  describes parallelization method for the coupled Maxwell - TDDFT calculations.

* Ref. :cite:`Hirokawa2016`
  describes computational aspects of electron dynamics calculations for periodic systems in many-core processors:

   
.. bibliography:: reference.bib
  :cited:
  :style: unsrt


