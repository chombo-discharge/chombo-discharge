.. PlasmaC documentation master file, created by
   sphinx-quickstart on Thu Jul 26 21:03:53 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**This is an alpha release of PlasmaC. Development is still in progress, and various bugs may or may not be present.**
   
Welcome to PlasmaC's user documentation!
========================================
	      
PlasmaC is a scalable computer code for two- and three-dimensional simulations of low-temperature fluid plasmas in complex geometries. PlasmaC features

* Multiphase Poisson and Helmholtz solvers
* Time-dependent multiphase heat solvers
* Diffusive radiative transport
* Convection-diffusion-reaction solvers

The solvers can be run on their own, or they can be coupled through one of our physics interfaces (for which PlasmaC was originally designed)

PlasmaC is built on top of `Chombo <https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations>`_, and therefore additionally features

* Cut-cell representation of multi-material geometries
* Patch based adaptive mesh refinement
* Scalable in parallel

Take a look at the :ref:`Chap:Gallery` for some examples. 

This documentation is the user documentation for usage and examples on how to use PlasmaC. There is also a separate :doxy:`Doxygen API <index>` documentation of the source code. 

Table of contents
=================

.. toctree::
   :maxdepth: 3

   Introduction
   Interface
   NewSimulations
   RunningSimulations
   ImportantClasses
   FAQ
   Bugs
   Todo
   Visualization
   Gallery

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _Chap:contact:

Contact us
==========

Send us an email at robert.marskar@sintef.no
