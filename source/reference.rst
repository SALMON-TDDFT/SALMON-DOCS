..
  We use sphinxcontrib-bibtex package for citing papers
  https://sphinxcontrib-bibtex.readthedocs.io/en/latest/index.html


.. _reference:

Reference
=================

Suggested Citations
--------------------
If you publish a paper in which SALMON makes an important contribution, please cite the SALMON code paper,
Ref. :cite:`Noda2019`
published in Computer Physics Communications.

We also suggest you to the following papers depending on your usage of SALMON.

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
