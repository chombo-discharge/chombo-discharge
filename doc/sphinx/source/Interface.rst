.. _Chap:Interface:

Understanding PlasmaC
=====================

PlasmaC is a comparatively flexible framework for plasma simulations. In fact, there are many abstractions in place that ensure that the code covers a broad range of physical phenomena, general geometries, multiple time stepping schemes, and general plasma-kinetic couplings. The mathematical form of this framework

.. math::

   \nabla\cdot\left(\epsilon_r\nabla\cdot\Phi\right) = -\frac{\rho}{\epsilon_0},

   \frac{\partial\sigma}{\partial t} = F_\sigma,

   \kappa\Psi - \nabla\cdot\left(\frac{1}{3\kappa}\nabla\Psi\right) = \frac{\eta}{c},

   \frac{\partial n}{\partial t} + \nabla\cdot\left(\mathbf{v} n - D\nabla n\right) = S.

This must be supported by additional boundary conditions on electrodes and insulating surfaces (on which the surface charge :math:`\sigma` lives). The number of advected species and radiative transport equations is arbitrary; the user can select the coupling through our interfaces: He can even turn off radiative transport altogether. 


PlasmaC works by embedding the equations above into an abstract C++ framework that the user must implement (or reuse existing pieces of) and compile into a mini-application. For most users, this will mostly include implementing a new geometry, a new plasma-kinetic scheme, or new functions for deciding when to coarsen or refine a certain spatial region. Internally, there are six such classes (described below), but most people will find it sufficient to modify at most three of these. 

.. toctree::
  :maxdepth: 1

  PlasmaCModules
  MiniApplications
  PythonInterface
  CodeStructure
