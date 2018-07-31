.. _Chap:ImportantClasses:

Class API
=========

Here, we discuss the base classes that makes up the foundation of PlasmaC.

PlasmaC uses a division-of-labor between physical and numerical modules. For a brief introduction to these, see :ref:`Chap:Interface`. The full, internal workings of PlasmaC are too complex to provide in detail, but we attempt to provide an overview here, as well as providing a summary of the input variables for the most important base classes.

Note that many of the classes below are abstract, and require implementation. For implementation-specific versions of these classes, please see the :ref:`Chap:Solvers` chapter. 

If you want to view the code in full, please see the :doxy:`Doxygen API <index>`.

Here are the base modules for PlasmaC, note that :ref:`Chap:plasma_kinetics`, :ref:`Chap:time_stepper`, :ref:`Chap:computational_geometry`, and :ref:`Chap:cell_tagger` are abstract, and require top-level implementation. 

.. toctree::
  :maxdepth: 1

  PlasmaEngine
  PlasmaKinetics
  AmrMesh
  PhysicalDomain
  TimeStepper
  ComputationalGeometry
  CellTagger
  GeoCoarsener
  Species
  Photon
