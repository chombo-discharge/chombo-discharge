.. _Chap:Interface:

Understanding PlasmaC
=====================

PlasmaC is a comparatively flexible framework for plasma simulations. In fact, there are many abstractions in place that ensure that the code covers a broad range of physical phenomena, general geometries, multiple time stepping schemes, and general plasma-kinetic couplings.

Broadly speaking, the user selects or implements a set of C++ abstraction classes that specify how the mini-app is built. For most users, this will mostly include implementing a new geometry, a new plasma-kinetic scheme, or new functions for deciding when to coarsen or refine a certain spatial region. Internally, there are six such classes (described below), but most people will find it sufficient to modify at most three of these.

.. toctree::
  :maxdepth: 1

  PlasmaCModules
  MiniApplications
  PythonInterface
  CodeStructure
